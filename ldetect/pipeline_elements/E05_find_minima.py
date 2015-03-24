import sys
import bisect
import ldetect.baselib.flat_file_consts as cnst
import ldetect.baselib.flat_file as flat
import ldetect.baselib.filters as filt
import ldetect.baselib.read_data as rd
import ldetect.baselib.binary_search as binsrch

class FlexibleBoundedAccessor:
	def __init__(self, data, f, min_ind, max_ind, invert):
		self.data = data
		self.f = f
		self.min_ind = min_ind
		self.max_ind = max_ind
		if invert:
			self.accessor = self.inverted_accessor
		else:
			self.accessor = self.regular_accessor

	def regular_accessor(self, i):
		if i >= self.min_ind and i<= self.max_ind:
			return self.f(self.data, i)
		else:
			raise Exception('Error: index out of range')

	def inverted_accessor(self, i):
		if i >= self.min_ind and i<= self.max_ind:
			return self.f(self.data, self.max_ind-i)
		else:
			raise Exception('Error: index out of range')

	def __getitem__(self, i):
		return self.accessor(i)

	def __len__(self):
		return self.max_ind-self.min_ind+1

# Moved to external binary search library:
# def find_le(a, x):
# 	'Find rightmost value less than or equal to x'
# 	i = bisect.bisect_right(a, x)
# 	if i:
# 		return i-1, a[i-1]
# 	raise ValueError

# def find_ge(a, x):
#     'Find leftmost item greater than or equal to x'
#     i = bisect.bisect_left(a, x)
#     if i != len(a):
#         return i, a[i]
#     raise ValueError

def find_end(data, f, x, val, max_srch_val=float('inf')):
	if x<=0:
		raise Exception('Error: x<=0')

	if x>=max_srch_val:
		raise Exception('Error: Max search value exceeded')

	if f(data, x) < val:
		return x
	else:
		return find_end(data, f, x*2, val)

def trackback(wrapper, srch_val, start_search, delta_coarse, step_coarse, step_fine=1):
	found_more = True # just to enter loop
	flat.print_log_msg('Starting coarse search')
	while found_more: # whenever we find more, search continues looking as far as delta width from newly found location
		found_more = False
		for i in range(start_search+step_coarse, start_search+delta_coarse, step_coarse):
			if wrapper[i] == srch_val:
				found_more = True
				start_search = i
				break

	# By default, step_fine = 1, therefore fine-grained search will happen. It can be disabled by setting step_fine to 0.
	if step_fine > 0:
		if step_fine > step_coarse:
			raise Exception('Error: step_fine is greater than step_coarse')
		delta_fine = step_coarse

		found_more = True # just to enter loop
		flat.print_log_msg('Starting fine search')
		while found_more: # whenever we find more, search continues looking as far as delta width from there
			found_more = False
			for i in range(start_search+step_fine, start_search+delta_fine, step_fine):
				if wrapper[i] == srch_val:
					found_more = True
					start_search = i		
					break

	return start_search

def custom_binary_search_with_trackback(np_init_array, f, srch_val, trackback_delta=200, trackback_step=20, init_search_location=1000):
	flat.print_log_msg('Starting custom_binary_search_with_trackback')

	# One-sided binary (i.e., exponential) search first
	end_v = find_end(np_init_array, f, init_search_location, srch_val)
	print('end_v: ', end_v)
	wrapper = FlexibleBoundedAccessor(np_init_array, f, 0, end_v, True)
	# Search with deferred detection of equality
	found_width_raw = binsrch.find_le_ind(wrapper, srch_val)
	found_width = end_v - found_width_raw
	print('found_width: ', found_width)

	# Find any remaining "noisy" minima
	found_width_trackback_raw = trackback(wrapper, srch_val, found_width_raw, trackback_delta, trackback_step)
	found_width_trackback = end_v - found_width_trackback_raw

	# Final result
	found_width = found_width_trackback
	print('found_width final: ', found_width)

	return found_width

def main():
	'''
	Run with at least:
	1 argument, defining the number of minima to be found

	optional last argument is file name of data to read into memory (if it hasn't alreay been read)
	'''

	if sys.argv[1] == 'help':
		print(main.__doc__)
		return

	srch_val = int(sys.argv[1])

	global np_init_array

	try:
		np_init_array		
	except NameError:
		flat.print_log_msg('Reading data')
		np_init_array, np_init_array_x = rd.read_data(sys.argv[len(sys.argv)-1])

	found_width = custom_binary_search_with_trackback(np_init_array, filt.apply_filter_get_minima, srch_val, trackback_delta=200, trackback_step=20, init_search_location=1000)
	print('found_width: ', found_width)
	flat.print_log_msg('Done')

if __name__ == '__main__':
	main()
