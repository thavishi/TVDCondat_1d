import numpy as np


def TVDCondat_1d(y, lamda):
	N = y.size;
	if N <= 1: 
		x = y 
		return x
	x             = np.zeros(len(y))# y can be a row or column vector.
	indstart_low  = np.zeros((N),dtype = int) # starting indices of constant segments of the lower approximation x_low
	indstart_up   = np.zeros((N), dtype = int) #starting indices of constant segments of the upper approximation x_up
	j_low         = 0 # index to count the segments of x_low
	j_up          = 0 # same for x_up
	jseg          = 0 # segment number of the current part under construction
	indjseg       = 0 # starting index of the current part

	#we have indjseg = indstart_low(jseg) = indstart_up(jseg)
	indstart_low[0]  = 0 # starting index of the j_low-th  segment of x_low
	indstart_up[0]   = 0 # same for x_up
	x_low_first      = y[0] - lamda   #value of the first segment of the part of x_low under construction
	x_up_first       = y[0] + lamda  # same for x_up
	x_low_curr       = x_low_first # value of the last segment of the part of x_low under construction
	x_up_curr        = x_up_first # same for x_up
	# the constant value of x_low over the j-th segment is stored in x(indstart_low(j_low)), except for j=jseg, where the value is x_low_first. Same for x_up. Indeed, the parts of  x_low and x_up under construction have distinct jump locations, but same starting index jseg.


	for i in range(1,N):
		# print(i)
		# print(x)
		if y[i] >= x_low_curr:
			if y[i] <= x_up_curr:
				#print(x_up_curr)
				#print(i)
				#print(j_up)
				#print(indstart_up)
				x_up_curr       = x_up_curr +(y[i] - x_up_curr)/(i - indstart_up[j_up]+1)
				#print(x_up_curr)
				x[indjseg]      = x_up_first
				while (j_up > jseg) and  (x_up_curr <= x[indstart_up[j_up-1]]):
					j_up      = j_up - 1
					x_up_curr = x[indstart_up[j_up]] + (x_up_curr - x[indstart_up[j_up]])*((i - indstart_up[j_up+1]+1)/(i-indstart_up[j_up]+1))

				if j_up == jseg:  # a jump in x downwards is possible       	
					# the fusion of x_low has not been done yet, but this is OK.
					# #print(x_up_curr)
					# #print(x_low_first)
					while (x_up_curr <= x_low_first) and (jseg < j_low):
			        	# the second test should always be true if the first one
			        	# is true and lamda>0, but this is a numerical safeguard.
			        	# And it is necessary if lamda=0.
			        	# validation of segments of x_low in x
						jseg                              = jseg + 1
						x[indjseg : indstart_low[jseg]]   = x_low_first
						x_up_curr                         = x_up_curr + (x_up_curr - x_low_first)*((indstart_low[jseg]-indjseg)/(i-indstart_low[jseg]+1))
						indjseg                           = indstart_low[jseg]
						x_low_first                       = x[indjseg]

					x_up_first        = x_up_curr
					j_up              = jseg
					indstart_up[jseg] = indjseg
				else: 
					x[indstart_up[j_up]]  = x_up_curr
			else: 
				# we start a new segment in x_up
				j_up               = j_up + 1
				indstart_up[j_up]  = i
				x[i]               = y[i]
				x_up_curr          = x[i]

			# fusion of x_low to keep it nonincreasing
			x_low_curr    = x_low_curr + (y[i] - x_low_curr)/(i - indstart_low[j_low] + 1)      
			x[indjseg]    = x_low_first

			while (j_low > jseg) and (x_low_curr >= x[indstart_low[j_low-1]]):
				j_low       = j_low - 1
				x_low_curr  = x[indstart_low[j_low]] + (x_low_curr - x[indstart_low[j_low]])*(( i - indstart_low[j_low+1]+1)/(i - indstart_low[j_low]+1 ))
			if j_low == jseg:  # a jump in x upwards is possible  
				while (x_low_curr >= x_up_first) and (jseg < j_up):
					# validation of segments of x_up in x
					jseg                           = jseg + 1
					x[indjseg:indstart_up[jseg]]   = x_up_first
					x_low_curr                     = x_low_curr + (x_low_curr - x_up_first)*((indstart_up[jseg]-indjseg)/(i-indstart_up[jseg]+1))
					indjseg                        = indstart_up[jseg]
					x_up_first                     = x[indjseg]
	       		
				x_low_first        = x_low_curr
				j_low              = jseg
				indstart_low[jseg] = indjseg
				if indjseg == i: # this part is not mandatory, it is a kind
					# of reset to increase numerical robustness.
					# If we are here, this is just after a jump upwards has
					# been validated. We have x_up_first=y(i).
					x_low_first = x_up_first - 2*lamda

			else: 
				x[indstart_low[j_low]] = x_low_curr
		else:
	    	# we start a new segment in x_low
			j_low                = j_low + 1
			indstart_low[j_low]  = i
			x[i]                 = y[i]
			x_low_curr           = x[i]
			# fusion of x_up to keep it nondecreasing
			x_up_curr            = x_up_curr + (y[i] - x_up_curr)/(i- indstart_up[j_up]+1)      
			x[indjseg]           = x_up_first

			while (j_up > jseg) and (x_up_curr <= x[indstart_up[j_up-1]]):
				j_up       = j_up - 1
				x_up_curr  = x[indstart_up[j_up]] + (x_up_curr - x[indstart_up[j_up]])*((i-indstart_up[j_up+1]+1)/(i-indstart_up[j_up]+1))

			if j_up == jseg:  # a jump in x downwards is possible 
				while (x_up_curr <= x_low_first) and (jseg < j_low): 
	        		# validation of segments of x_low in x
					jseg                             = jseg + 1
					x[indjseg:indstart_low[jseg]]    = x_low_first
					x_up_curr                        = x_up_curr + (x_up_curr - x_low_first)*((indstart_low[jseg] - indjseg)/(i - indstart_low[jseg] + 1))
					indjseg                          = indstart_low[jseg]
					x_low_first                      = x[indjseg]
	   			
				x_up_first        = x_up_curr
				j_up              = jseg
				indstart_up[jseg] = indjseg
				if indjseg == i: # this part is not mandatory, it is a kind
	       			# of reset to increase numerical robustness.
					x_up_first = x_low_first + 2*lamda
			else: 
				x[indstart_up[j_up]] = x_up_curr

	i = N-1
	if y[i] + lamda <= x_low_curr: 
		# the segments of x_low are validated
		while jseg < j_low:
			jseg                             = jseg + 1
			x[indjseg:indstart_low[jseg]]    = x_low_first
			indjseg                          = indstart_low[jseg]
			x_low_first                      = x[indjseg]

		x[indjseg:i] = x_low_first
		x[i]           = y[i] + lamda
	elif y[i]-lamda >= x_up_curr: 
		# the segments of x_up are validated
		while jseg < j_up:
			jseg                            = jseg + 1
			x[indjseg:indstart_up[jseg]]    = x_up_first
			indjseg                         = indstart_up[jseg]
			x_up_first                      = x[indjseg]
	 	
		x[indjseg:i] = x_up_first
		x[i]           = y[i] - lamda
	else:
		# fusion of x_low to keep it nonincreasing
		x_low_curr = x_low_curr + (y[i] + lamda - x_low_curr)/(i - indstart_low[j_low] + 1)      
		x[indjseg] = x_low_first
		while (j_low > jseg) and (x_low_curr >= x[indstart_low[j_low-1]]):
			j_low      = j_low - 1
			x_low_curr = x[indstart_low[j_low]] + (x_low_curr - x[indstart_low[j_low]])*((i-indstart_low[j_low+1]+1)/(i-indstart_low[j_low]+1))

		if j_low == jseg: # the segments of x_up must be validated
			if x_up_first >= x_low_curr: # same unique segment of x_low and x_up
				x[indjseg:i+1] = x_low_curr
			else:
				# fusion of x_up to keep it nondereasing
				x_up_curr   = x_up_curr + (y[i] - lamda - x_up_curr)/(i - indstart_up[j_up] + 1)      
				x[indjseg]  = x_up_first
				while (j_up > jseg) and (x_up_curr <= x[indstart_up[j_up-1]]):
					j_up      = j_up - 1
					x_up_curr = x[indstart_up[j_up]] + (x_up_curr - x[indstart_up[j_up]])*((i - indstart_up[j_up+1]+1)/(i - indstart_up[j_up] + 1))

					x[indstart_up[j_up]:i+1] = x_up_curr
				while jseg < j_up:  # the segments of x_up are validated
					jseg                           = jseg + 1
					x[indjseg:indstart_up[jseg]]   = x_up_first
					indjseg                        = indstart_up[jseg]
					x_up_first                     = x[indjseg]
		else: 	# the segments of x_low must be validated
			x[indstart_low[j_low]:i+1] = x_low_curr
			while jseg < j_low:
				jseg                            = jseg + 1
				x[indjseg:indstart_low[jseg]]   = x_low_first
				indjseg                         = indstart_low[jseg]
				x_low_first                     = x[indjseg]

	return x