class Collector(object):

	def __init__(self, base_name=None):
		self.base_name = base_name
		self.tool_name = ''
		self.version = 0
		return

	def parse(self, file_path, descriptor=None):
		metrics_file = open(file_path, 'r')
		metrics_dict = {}

		for line in metrics_file:
			line = line.lower()

			try:
				key, value = line.split('|')
				key = key.lstrip().strip()
				key = key.replace('%', 'percent')
				key = key.replace(' ', '_')
				key = key.replace(':', '')
				key = key.replace('(', '')
				key = key.replace(')', '')
				key = key.replace('_of', '')
				key = key.strip(',')[0]

				value = value.lstrip().strip()
				value = value.replace('%', '')
				metrics_dict[key] = value

			except:
				pass

		metrics_dict.pop('started_mapping_on', None)
		metrics_dict.pop('started_job_on', None)
		metrics_dict.pop('finished_on', None)

		return metrics_dict

	def record_metrics(self):
		pass

if __name__ == '__main__':

	import sys
	
	collector = Collector()
	for key in collector.parse(sys.argv[0]):
		print key
