from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import glob
import compare_two_zlists as cv
import math
from scipy.stats import hypergeom
from decimal import Decimal
from math import log 

def make_venn_matrix(filename_list):
	fig1 = plt.figure(1)
	fig1.suptitle('Differentially Expressed Genes Overlap', fontsize=24)

	subplot_counter = 1

	print len(filename_list)*len(filename_list)
	for zlist1 in filename_list:
		for zlist2 in filename_list:

			plt.subplot(len(filename_list), len(filename_list), subplot_counter)
	
		
			offset = math.ceil(float(subplot_counter)/float(len(filename_list)))
			position = int(subplot_counter - 1) % len(filename_list) + 1

			if zlist1 == zlist2:

				plt.subplot(len(filename_list), len(filename_list), subplot_counter)

				up_A, up_B, all_changing = cv.get_changing(zlist1)
				plt.text(0, 0, '{}'.format(zlist1))
				plt.text(0, .4, 'Higher in A: {}'.format(str(len(up_A))))
				plt.text(0, .2, 'Higher in B: {}'.format(str(len(up_B))))

				plt.axis('off')
				plt.plot()

				print 'working {}'.format(subplot_counter)
				subplot_counter+=1

			else:
				plt.subplot(len(filename_list), len(filename_list), subplot_counter)
				
				(venn_values, all_union) = cv.compare_two(zlist1, zlist2)

				color1 = ''
				color2 = ''

				if position > offset:

					color1 = 'MediumVioletRed'
					color2 = 'OrangeRed'
					union = venn_values['up_A']['union']
					in_common = venn_values['up_A']['common']
					unique_1 = venn_values['up_A']['up_1']
					unique_2 = venn_values['up_A']['up_2']

				if position < offset:
					
					color1 = 'LimeGreen'
					color2 = 'DodgerBlue'
					union = venn_values['up_B']['union']
					in_common = venn_values['up_B']['common']
					unique_1 = venn_values['up_B']['up_1']
					unique_2 = venn_values['up_B']['up_2']

				total_genes = len(all_union)
				total_1 = unique_1 + in_common
				total_2 = unique_2 + in_common

				try:
					log_prob = Decimal(log(hypergeom.sf(in_common, total_genes, total_1, total_2)))

				except:
					log_prob = '-inf'

				plt.plot(cv.draw_venn(union, in_common, unique_1, unique_2, color1, color2))

				if log_prob != '-inf':
					plt.annotate('log p-value: %2.3f'%log_prob, xy=(0,0), xycoords='axes fraction')
				else:
					plt.annotate('log p-value: -inf', xy=(0,0), xycoords='axes fraction')
				print 'working {}'.format(subplot_counter)
				print str(total_genes)
				print str(log_prob)
				subplot_counter+=1
				
	plt.show()
	return

if __name__ == '__main__':
	
	venn_filelist = glob.glob('*zlist')
	venn_filelist.sort()
	make_venn_matrix(venn_filelist)

