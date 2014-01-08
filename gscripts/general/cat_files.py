import argparse
import os
from glob import glob
import shutil
from collections import defaultdict
from copy import copy

def inspectFolder(folder, extensions, recurse=False, cat=False, cat_depth=1):
    ext_list = []
   
    if isinstance(extensions, basestring):
        ext_list.append(extensions)
    
    else:
        ext_list = extensions
  
    if os.path.isdir(folder):
        os.chdir(folder)
        
        for ext in ext_list:
            
            file_list = glob('*{}'.format(ext))
            
            if file_list:
                print
                print folder
                print 'extension: {}'.format(ext)
                if cat:
                    catFiles(file_list, ext, cat_depth)
                
                
            
        for file in glob('*'):
                if recurse:
                    if os.path.isdir(file):
                        inspectFolder(file, extensions, recurse, cat, cat_depth)
                        
        os.chdir('..')
            
    else:
        print '{} is not a directory.'.format(str(folder))
        
    return


def catFiles(file_list, ext, combine_depth=1):
    
    # ASSUMPTION: file name VALUES go from more general to more specific
    # ex. SAMPLE_INDEX_LANE_MATE_FILENO.fastq.gz
    
    
    if isinstance(file_list, basestring):
        return

    filenames_tree = defaultdict(defaultdict)
    
    for file in file_list:
        complete_file = file
        file = file.split('.')[0]
        
        file_values = file.split('_')
        
        current_parent = None
        current_file_string = ''
        for step in file_values:
            
            if current_file_string:
                current_file_string += '_{}'.format(step)
                
            else:
                current_file_string = step
            if step == file_values[0]:
                
                if not step in filenames_tree:
                    filenames_tree[step] = defaultdict(defaultdict) 
                    filenames_tree[step]['file_string'] = current_file_string
                    
                current_parent = filenames_tree[step]
                    
            else:
                
                if step != file_values[-1]:
                    
                    if not step in current_parent:
                        current_parent[step] = defaultdict(defaultdict)
                        
                    current_parent[step]['file_string'] = current_file_string
                    current_parent = current_parent[step]
                    
                    
                else:
                    current_parent[step] = complete_file
                    
        
    for item in gatherTrees(filenames_tree, combine_depth, 0, []):
        
        for key in item:
            if not key == 'file_string':
                
                destination_file = '_'.join([str(item['file_string']), str(key)])
                destination_file += ext
                
                if destination_file in file_list:
                    raise Exception('Overwrite of original file {}'.format(destination_file))
                leaf_list = recurseTree(item[key], [])
                leaf_list.sort()
                
                print
                print 'cat to:'
                print destination_file
                print
                print 'the files:'
                for some_file in leaf_list:
			print some_file

		if not os.path.exists(destination_file):
			destination = open(destination_file, 'wb')

			for leaf in leaf_list:
				shutil.copyfileobj(open(leaf, 'rb'), destination)

			destination.close()


def gatherTrees(tree, target_depth=1, depth=0, tree_list=[]):
    depth += 1
    
    if target_depth == depth:
        tree_list.append(tree)
    
    if isinstance(tree, defaultdict):
        for key in tree:
            if isinstance(tree[key], defaultdict):
                gatherTrees(tree[key], target_depth, depth, tree_list)
            
    else:
        return tree
    
    return tree_list


def recurseTree(tree, leaf_list=[]):
    
    if isinstance(tree, defaultdict):
        
        for key in tree:
            if not key == 'file_string':
                recurseTree(tree[key], leaf_list)
             
    elif isinstance(tree, basestring):
        #print tree
        leaf_list.append(tree)
        
    else:
        return tree
        
    return leaf_list
        


#inspectFolder('/nas3/yeolab/SolexaData/20130430_Single_cell_RNASeq/', '.fastq.gz', True, True, 5)


if __name__ == '__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--recurse", action="store_true")
    parser.add_argument("-c", "--cat", action="store_true")
    parser.add_argument("-d", "--directory", dest="directory")
    parser.add_argument("-e", "--extension", dest="extension")
    parser.add_argument("-p", "--depth", type=int, dest="depth")
    args = parser.parse_args()

    inspectFolder(args.directory, args.extension, args.recurse, args.cat, args.depth)



