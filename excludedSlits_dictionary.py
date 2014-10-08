#!/usr/bin/env python
import numpy
#coding: utf8 
#Dictionary of the slits to exclude from the galaxies in the sample



#Points excluded from the dataset (indices)

					#Galaxy		#index	#ID slit 	#RA				#Dec 					#Reason of exclusion
dicExcludedSlits = {#'NGC1023': [[128, 	'hst_260',  20.549987411508, 8.370000762947999],	# too low CaT index near the galaxy centre
					#			], 
					'NGC5846': [[3,	'add_GC_0126361',  -10.8600080109, -33.1500015259],	# Point within 15 arcsec from NGC5846a
					[48,	'add_GC_0125285',  -12.8099954224, -43.8499984741],	# Point within 15 arcsec from NGC5846a
					[69,	'add_GC_0125763',  3.24, -28.5499992371],	# Point within 15 arcsec from NGC5846a
					[70,	'add_GC_0126397', -2.90999771127, -27.4500007629],	# Point within 15 arcsec from NGC5846a
					[72, 'add_GC_0126397', -2.90999771127, -27.4500007629],	# Point within 15 arcsec from NGC5846a
					[92, 'add_GC_0125588', 4.43999885562, -29.1500015259],	# Point within 15 arcsec from NGC5846a
					[94, 'add_GC_0126360', -9.81001258852, -32.3499984741],					# Point within 15 arcsec from NGC5846a
					[78, 'scam_0138168',	-57.65999198904, 46.5500001906],	# Extreme point 
						],
					'NGC3607': [[161, 'SKiMS_2', -133.335, -117.42000038136]	# Too far
					]
					}
