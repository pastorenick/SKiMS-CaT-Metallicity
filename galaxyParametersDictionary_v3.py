#!/usr/bin/env python
import numpy
#coding: utf8 
#Dictionaries for the galaxies in the sample
# v3 -> Match with Jacob's table, except some ellipticities

#From Cappellari+11a, except N720, N1400, N1407, N3115 derived using the same methodology
# NGC4449 -> Rich+2012
Reff = {"NGC720": 33.9, "NGC821":39.81, "NGC1023":47.86, "NGC1400": 29.3,
        "NGC1407":63.4, "NGC2768":63.10, "NGC2974":38.02, "NGC3115":32.1,
        "NGC3377":35.48, "NGC3607":38.9045, "NGC3608":29.51, "NGC4111":12.02, "NGC4278":31.62, 
        "NGC4365":52.48, "NGC4342":6.60693, "NGC4374":52.48, "NGC4473":26.92, "NGC4486":81.28, 
        "NGC4449": 99.33,
        "NGC4494": 48.98, "NGC4526":44.67, "NGC4564":20.4174, "NGC4649":66.07, "NGC4697":61.66, 
        "NGC5846":58.88, "NGC5866":36.31, 
        "NGC5907":20,	#FAKE...
        "NGC7457":36.31}

Sigma0_Unknown = {"NGC720": 241., "NGC821":189., "NGC1023": 182.,"NGC1400":254. ,
        "NGC1407":279., "NGC2768":181., "NGC2974":233., "NGC3115":230.,
        "NGC3377":138., "NGC3607":224., "NGC4111":151., "NGC4278":231., 
        "NGC4365":265., "NGC4374":283., "NGC4473":180., "NGC4486":298., 
        "NGC4494":145., "NGC4526":222.,"NGC4649":341., "NGC4697":171., 
        "NGC5846":238., "NGC5866":172., "NGC7457":78.}

#The adopted sigmas, all from HyperLeda
Sigma0_HL = {"NGC720": 241., "NGC821":200., "NGC1023": 204.,"NGC1400":252. ,
        "NGC1407":271., "NGC2768":181., "NGC2974":238., "NGC3115":267.,
        "NGC3377":139., "NGC3607":224., "NGC3608":192.,
        "NGC4111":149., "NGC4278":237., "NGC4342":252., 
        "NGC4365":256., "NGC4374":283., "NGC4473":179., "NGC4486":334., 
        "NGC4449":18., 
        "NGC4494":150., "NGC4526":251., "NGC4564":157.,
        "NGC4649":335., "NGC4697":171., 
        "NGC5846":239., "NGC5866":159., "NGC5907":121., "NGC7457":69.}

CentreCoordinates = {"NGC720": ['01:53:00.498','-13:44:19.15'], "NGC821":['02:08:21.14','+10:59:41.7'], 
  "NGC1023":['02:40:24.010','+39:03:47.83'], "NGC1400":['03:39:30.836','-18:41:17.05'],
  "NGC1407":['03:40:11.86','-18:34:48.4'], "NGC2768":['09:11:37.50','+60:02:14.0'], 
  "NGC2974":['09:42:33.28','-03:41:56.9'], "NGC3115":['10:05:13.98','-07:43:06.9'],
  "NGC3377":['10:47:42.331','+13:59:09.30'], "NGC3607":['11:16:54.639','+18:03:06.32'], 
  "NGC3608":['11:16:58.951','+18:08:55.26'], 
  "NGC4111":['12:07:03.132','+43:03:56.59'], "NGC4278":['12:20:06.8242','+29:16:50.722'],
  "NGC4342":['12:23:39.000','+07:03:14.38'],
  "NGC4365":['12:24:28.284','+07:19:03.62'], "NGC4374":['12:25:03.7433','+12:53:13.139'], 
  "NGC4473":['12:29:48.870','+13:25:45.69'], "NGC4486":['12:30:49.423','+12:23:28.04'], 
  "NGC4449":['12:28:11.103','+44:05:37.07'], 
  "NGC4494":['12:31:24.104','+25:46:30.91'], "NGC4526":['12:34:03.085','+07:41:58.29'], 
  "NGC4564":['12:36:26.983','+11:26:21.42'], 
  "NGC4649":['12:43:39.975','+11:33:09.74'], "NGC4697":['12:48:35.878','-05:48:02.67'], 
  "NGC5846":['15:06:29.284','+01:36:20.25'], "NGC5866":['15:06:29.4988','+55:45:47.568'], 
  "NGC5907":['15:15:53.770','+56:19:43.58'], 
  "NGC7457":['23:00:59.934','+30:08:41.79']}


# From Krajnovic+11, except NGC720 -> Cappellari+07, N1400 and N1407 -> Spolaor08, 
# NGC4449 -> RC3
PA0 = {"NGC720": 142.3, "NGC821":31.2,  "NGC1023":83.3,"NGC1400":36.1,
        "NGC1407":58.3, "NGC2768":91.6, "NGC2974":44.2, "NGC3115":43.5,
        "NGC3377":46.3, "NGC3607":124.8, "NGC3608":82.0, "NGC4111":150.3, "NGC4278":39.5, 
        "NGC4365":40.9, "NGC4342":163.6, "NGC4374":128.8, "NGC4473":92.2, 
        "NGC4449":44.0, 
        "NGC4486":151.3, 
        "NGC4494":176.3, "NGC4526":113.7, "NGC4564":48.5, "NGC4649":91.3, "NGC4697":67.2, 
        "NGC5846":53.3, "NGC5866":125.0, 
        "NGC5907":155., #From 2Mass
        "NGC7457":124.8}

# From Emsellem+11, except NGC720 -> Cappellari+07 -> Cappellari+07, N1400 and N1407 -> Spolaor08
# NGC3115 Capaccioli+87, NGC4449 -> RC3
b_a = {"NGC720": 0.57, "NGC821":0.65, "NGC1023":0.63,"NGC1400":0.89,
        "NGC1407":0.95, "NGC2768":0.53, "NGC2974":0.59, "NGC3115":0.51,
        "NGC3377":0.50, "NGC3607":0.87, "NGC3608":0.80, "NGC4111":0.42, "NGC4278":0.90, 
        "NGC4365":0.78, "NGC4342":0.42, "NGC4374":0.85, 
        "NGC4449":0.71, "NGC4473":0.58, "NGC4486":0.96, 
        "NGC4494":0.83, "NGC4526":0.64, "NGC4564":0.47, "NGC4649":0.84, "NGC4697":0.55, 
        "NGC5846":0.94, "NGC5866":0.43, 
        "NGC5907":0.16, #From 2Mass
        "NGC7457":0.53}

#From NED
vel0 = {"NGC720": 1745., "NGC821":1718., "NGC1023": 602.,"NGC1400":558.,
        "NGC1407":1779., "NGC2768":1353., "NGC2974":1919., "NGC3115":663.,
        "NGC3377":690., "NGC3607":960., "NGC3608":1226., "NGC4111":792., "NGC4278":620., 
        "NGC4365":1243., "NGC4374":1017., "NGC4473":2244., "NGC4486":1284., 
        "NGC4494":1342., "NGC4526":617., 
        "NGC4564":1142., "NGC4649":1110., "NGC4697":1252., 
        "NGC5846":1714., "NGC5866":672., "NGC5907":667., "NGC7457":844.}


#From Jacob's
Morph = {"NGC720": 'E', "NGC821":'E', "NGC1023":'S0', "NGC1400": 'S0',
        "NGC1407":'E', "NGC2768":'E', "NGC2974":'E', "NGC3115":'S0',
        "NGC3377":'E', "NGC3607":'S0', "NGC4111":'S0', "NGC4278":'E', 
        "NGC4365":'E', "NGC4374":'E', "NGC4473":'E', "NGC4486":'E', 
        "NGC4494": 'E', "NGC4526":'S0', "NGC4649":'E', "NGC4697":'E', 
        "NGC5846":'E', "NGC5866":'S0', "NGC7457":'S0'}

#All from ATLAS3d excluding N720, N1400, N1407 and N3115 (which are taken from Jacob's)
distMpc = {"NGC720": 26.9, "NGC821":23.4, "NGC1023":11.1, "NGC1400":26.8, 
        "NGC1407":26.8, "NGC2768":21.8, "NGC2974":20.9, "NGC3115":9.4,
        "NGC3377":10.9, "NGC3607":22.2, "NGC4111":14.6, "NGC4278":15.6, 
        "NGC4365":23.3, "NGC4374":18.5, "NGC4473":15.3, "NGC4486":17.7, 
        "NGC4494": 16.6, "NGC4526":16.4, "NGC4649":17.3, "NGC4697":11.4, 
        "NGC5846":24.2, "NGC5866":14.9, "NGC7457":12.9}

#All from ATLAS3d excluding N720, N1400, N1407 and N3115 (which are taken from Jacob's)
LK_adopted = {"NGC720":-24.6 , "NGC821":-23.99, "NGC1023":-24.01, "NGC1400":-24.3,
        "NGC1407":-25.4, "NGC2768":-24.71, "NGC2974":-23.62, "NGC3115":-24.0,
        "NGC3377":-22.76, "NGC3607":-24.74, "NGC4111":-23.27, "NGC4278":-23.80, 
        "NGC4365":-25.21, "NGC4374":-25.12, "NGC4473":-23.77, "NGC4486":-25.38, 
        "NGC4494": -24.11, "NGC4526":-24.62, "NGC4649":-25.46, "NGC4697":-23.93, 
        "NGC5846":-25.01, "NGC5866":-24.00, "NGC7457":-22.38}



LK_2mass = {"NGC720": [7.271, 0.024], "NGC821": [7.900,0.019], "NGC1023": [6.238,0.021],"NGC1400":[7.811,0.018],
        "NGC1407":[6.702, 0.025], "NGC2768": [6.997, 0.031], "NGC2974": [6.255, 0.005], "NGC3115": [5.883, 0.017],
        "NGC3377": [7.441, 0.029], "NGC3607":[6.994, 0.023], "NGC4111": [7.553, 0.018], "NGC4278": [7.184, 0.011], 
        "NGC4365":[6.640,0.029], "NGC4374":[6.222,0.023], "NGC4473":[7.157,0.023], "NGC4486":[5.812,0.019], 
        "NGC4494":[6.998,0.024], "NGC4526":[6.473,0.020], "NGC4649":[5.739,0.021], "NGC4697":[6.367,0.027], 
        "NGC5846":[6.935,0.023], "NGC5866":[6.873,0.018], "NGC7457":[8.192,0.024]}


#From NED
HelioVel_OLD = {"NGC720": [1620.0,31.0], "NGC821": [1735.0,4.0], "NGC1023": [637.0, 4.0],"NGC1400":[530.0,45.],
        "NGC1407":[1779.,9.], "NGC2768": [1373.,5.], "NGC2974": [1919.,13.], "NGC3115": [663.0,4.],
        "NGC3377": [665.,2.], "NGC3607":[960.,20.], "NGC4111": [807.,10.], "NGC4278": [649.,5.], 
        "NGC4365":[1243.,6.], "NGC4374":[1060.,6.], "NGC4473":[2244.,2.], "NGC4486":[1307.,7.], 
        "NGC4494":[1344.,11.], "NGC4526":[448.,8.], "NGC4649":[1117.,6.], "NGC4697":[1241.,2.], 
        "NGC5846":[1630.,18.], "NGC5866":[672.,9.], "NGC7457":[812.,6.]}

#From Cappellari+11 except N720, N1400, N1407, N3115, N4449, which are from NED
HelioVel = {"NGC720": 1745., "NGC821": 1718., "NGC1023": 602.,"NGC1400":558.,
        "NGC1407":1779., "NGC2768": 1353., "NGC2974": 1887., "NGC3115": 663.,
        "NGC3377": 690., "NGC3607":942., "NGC3608":1226., "NGC4111":792., "NGC4278": 620., 
        "NGC4449": 207.,
        "NGC4365":1243., "NGC4374":1017., "NGC4473":2260., "NGC4486":1284., 
        "NGC4494":1342., "NGC4526":617., "NGC4564":1142., "NGC4649":1110., "NGC4697":1252., 
        "NGC5846":1712., "NGC5866":755., "NGC5907":667., "NGC7457":844.}



galLiterature_metallicity_values = {	#PAPER, x1, y1,[erry1+, erry1-], x2, y2, ..., Reff1, Reff2, ...
	"NGC720": 	[	['Trager et al. (2000a)', 1./8., 0.44, [0.15,0.15],	1./2., 1.13, [0.42, 0.42],	40.]
				,['Howell (2005)', 1./8., 0.55, [0.09, 0.09],	66.]
			],
	"NGC821": 	[	['Trager et al. (2000a)', 1./8., 0.22, [0.03, 0.03],	1./2.,	0.12, [0.05, 0.05],	36.]
				,['Howell (2005)', 1./8., 0.55, [0.09, 0.09],	66.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.48, [0.13, 0.17],		50.]
				,['Conroy et al. (2012)', 1./8., 0.14, [0.,0.],		31.3]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC821"], 0.11, [0.03, 0.03], 	Reff["NGC821"]]	#Within 4arcsec
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.46, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
			],
	"NGC1023": 	[	['Conroy et al. (2012)', 1./8., 0.16, [0.,0.],		48.7]
				,[r"Sil'chenko (2006)", 4./Reff['NGC1023'], 0.4, [0.,0.],		7./Reff['NGC1023'], -0.1, [0.,0.], Reff['NGC1023']]
			],
	"NGC1400": 	[	['Spolaor et al. (2008b)', 1./8., 0.25, [0.06, 0.06],	27.]
				,['Howell (2005)', 1./8., 0.31, [0.13, 0.13],	63.]
				,[r"Barr et al. (2007)", 1./8., 0.67, [0.04, 0.04],	29.32]
				,["Idiart et al. (2007)", 1./8., 0.02, [0., 0.],	63.]
			],
	"NGC1407": 	[	['Spolaor et al. (2008b)', 1./8., 0.29, [0.08, 0.08],	72.]
				,['Howell (2005)', 1./8., 0.56, [0.07, 0.07],	120.]
				,[r"Cenarro et al. (2007)", 1./8., 0.33, [0.03, 0.03],	36.4]
				,["Annibali et al. (2007)", 1./8., 0.03, [0.01, 0.01],	70.3]
			],
	"NGC2768": 	[	['Howell (2005)', 1./8., 0.14, [0.25, 0.25],	89.]
				,[r"Denicol\'{o} et al. (2005b) - Major", 1./8., 0.38, [0.29, 0.30],		64.]
				,[r"Denicol\'{o} et al. (2005b) - Minor", 1./8., 0.64, [0.06, 0.03],		64.]
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC2768'], 1.38, [0.,0.], 	Reff['NGC2768']]	#Within 2.5arcsec
				,['Conroy et al. (2012)', 1./8., 0.09, [0.,0.],		68.]
				,["Serra et al. (2008)", 1./16., 0.39, [0.04,0.04],		1./2., 0.27, [0.03,0.06],	27.]
				,[r"Sil'chenko (2006)", 4./Reff['NGC2768'], 0.1, [0.,0.],		7./Reff['NGC2768'], 0., [0.,0.], Reff['NGC2768']]
				,["Idiart et al. (2007)", 1./8., 0.18, [0., 0.],	89.]
			],
	"NGC2974": 	[	[r"Denicol\'{o} et al. (2005b)", 1./8., 0.15, [0.06, 0.07],		24.]
				,['Conroy et al. (2012)', 1./8., 0.13, [0.,0.],		28.3]
				,["Annibali et al. (2007)", 1./8., 0.02, [0.01, 0.01],	24.4]
			],
	"NGC3115": 	[	['Howell (2005)', 1./8., 0.65, [0.06, 0.06],	60.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC3115"], 0.13, [0.03, 0.03], 	Reff["NGC3115"]]	#Within 4arcsec
				,["Idiart et al. (2007)", 1./8., 0.36, [0., 0.],	60.]
			],
	"NGC3377": 	[	['Trager et al. (2000a)', 1.8, 0.19, [0.06, 0.06],		1./2.,	-0.12, [0.04, 0.04],	34.]
				,['Howell (2005)', 1./8., 0.30, [0.06, 0.06],	56.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.35, [0.05, 0.05],		34.]
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC3377'], 1.71, [0.,0.], 	Reff['NGC3377']]	#Within 2.5arcsec
				,['Conroy et al. (2012)', 1./8., 0.12, [0.,0.],		38.3]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC3377"], -0.01, [0.07, 0.07], 	Reff["NGC3377"]]	#Within 4arcsec
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.41, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
				,["Idiart et al. (2007)", 1./8., -0.02, [0., 0.],	56.]
				#,["Harris et al. (2007)", 99./Reff['NGC3377']., -0.6, [0.3., 0.9],	56.]#From stellar halo measures
			],
	"NGC3607": 	[	['Howell (2005)', 1./8., 0.27, [0.10, 0.10],	109.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.67, [0.12, 0.07],		43.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC3607"], 0.01, [0.06, 0.06], 	Reff["NGC3607"]]	#Within 4arcsec
				,['Proctor et al. (2002)', 2.12/Reff['NGC3607'], 0.49, [0.05,0.05],		Reff['NGC3607']]	#Central area 3.6x1.25 arcsec 
				,[r"Sil'chenko (2006)", 4./Reff['NGC3607'], 0.2, [0.,0.],		7./Reff['NGC3607'], 0.2, [0.,0.], Reff['NGC3607']]
				,["Idiart et al. (2007)", 1./8., 0.05, [0., 0.],	109.]
				,["Annibali et al. (2007)", 1./8., 0.05, [0.01, 0.01],	43.4]
			],
	"NGC4111": 	[	[r"Sil'chenko (2006)", 4./Reff['NGC4111'], 0.6, [0.,0.],		7./Reff['NGC4111'], 0., [0.,0.], Reff['NGC4111']]
			],
	"NGC4278": 	[	['Conroy et al. (2012)', 1./8., 0.05, [0.,0.],		30.6]
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC4278'], 2.33, [0.,0.], 	Reff['NGC4278']]	#Within 2.5arcsec
				,["Serra et al. (2008)", 1./16., 0.60, [0.07,0.11],		1./2., 0.59, [0.04,0.05],	26.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC4278"], 0.08, [0.03, 0.03], 	Reff["NGC4278"]]	#Within 4arcsec
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.30, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
			],
	"NGC4365": 	[	['Howell (2005)', 1./8., 0.59, [0.08, 0.08],	95.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.65, [0.10, 0.06],		50.]
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC4365'], 2.00, [0.,0.], 	Reff['NGC4365']]	#Within 2.5arcsec
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC4365"], 0.14, [0.03, 0.03], 	Reff["NGC4365"]]	#Within 4arcsec
				,['Proctor et al. (2002)', 2.12/Reff['NGC4365'], 0.41, [0.04,0.04],		Reff['NGC4365']]	#Central area 3.6x1.25 arcsec 
				,["Idiart et al. (2007)", 1./8., 0.36, [0., 0.],	95.]
			],
	"NGC4374": 	[	['Trager et al. (2000a)', 1.8, 0.12, [0.03, 0.03],		1./2.,	-0.01, [0.05, 0.05],	52.]
				,['Howell (2005)', 1./8., 0.24, [0.02, 0.02],	91.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.50, [0.11, 0.15],		51.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC4374"], 0.10, [0.03, 0.03], 	Reff["NGC4374"]]	#Within 4arcsec
				,["Gallagher et al. (2008)", 9.63/Reff['NGC4374'], 0.10, [0.06, 0.05],		Reff['NGC4374']]	#Within 9.63arcsec
				,['Proctor et al. (2002)', 2.12/Reff['NGC4374'], 0.28, [0.06,0.06],		Reff['NGC4374']]	#Central area 3.6x1.25 arcsec 
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.38, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
				,["Idiart et al. (2007)", 1./8., 0.22, [0., 0.],	91.]
				,["Annibali et al. (2007)", 1./8., 0.03, [0.01, 0.01],	50.9]
			],
	"NGC4473": 	[	['Howell (2005)', 1./8., 0.56, [0.08, 0.08],	42.]
				,['Conroy et al. (2012)', 1./8., 0.18, [0.,0.],		28.6]
				,["Idiart et al. (2007)", 1./8., 0.24, [0., 0.],	42.]
			],
	"NGC4494": 	[	[r"Denicol\'{o} et al. (2005b)", 1./8., 0.19, [0.04, 0.04],		49.]
			],
	"NGC4526":	[	["Gallagher et al. (2008)", 10.78/Reff['NGC4526'], 0.40, [0.06, 0.06],		Reff['NGC4526']]	#Within 10.78arcsec
				,['Proctor et al. (2002)', 2.12/Reff['NGC4526'], 0.49, [0.06,0.06],		Reff['NGC4526']]	#Central area 3.6x1.25 arcsec 
				,[r"Sil'chenko (2006)", 4./Reff['NGC4526'], 0.4, [0.,0.],		7./Reff['NGC4526'], 0.1, [0.,0.], Reff['NGC4526']]
			],
	"NGC4649": 	[	['Trager et al. (2000a)', 1.8, 0.27, [0.04, 0.04],		1./2.,	0.05, [0.04, 0.04],	74.]
				,['Howell (2005)', 1./8., 0.37, [0.03, 0.03],	123.]
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC4649'], 2.50, [0.,0.], 	Reff['NGC4649']]	#Within 2.5arcsec
				,["Gallagher et al. (2008)", 13.09/Reff['NGC4649'], -0.08, [0.05, 0.06],		Reff['NGC4649']]	#Within 13.09arcsec
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.58, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
			],
	"NGC4697": 	[	['Trager et al. (2000a)', 1.8, 0.06, [0.05, 0.05],		1./2.,	-0.29, [0.04, 0.04],	75.]
				,['Howell (2005)', 1./8., 0.19, [0.04, 0.04],	126.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC4697"], 0.11, [0.03, 0.03], 	Reff["NGC4697"]]	#Within 4arcsec
				,['Proctor et al. (2002)', 2.12/Reff['NGC4697'], 0.49, [0.05,0.05],		Reff['NGC4697']]	#Central area 3.6x1.25 arcsec 
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.33, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
				,["Idiart et al. (2007)", 1./8., -0.02, [0., 0.],	126.]
				,["Annibali et al. (2007)", 1./8., 0.02, [0.00, 0.00],	72.]
			],
	"NGC5846": 	[	['Trager et al. (2000a)', 1./8., 0.14, [0.05, 0.05],		1./2.,	-0.13, [0.04, 0.04],	83.]
				,['Howell (2005) - Major', 1./8., 0.44, [0.14, 0.14],	138.]
				,['Howell (2005) - Minor', 1./8., 0.25, [0.04, 0.04],	138.]
				,[r"Denicol\'{o} et al. (2005b)", 1./8., 0.34, [0.17, 0.18],		62.]
				,[r"S{\'a}nchez-Bl{\'a}zquez et al. (2006)", 4./Reff["NGC5846"], 0.10, [0.03, 0.03], 	Reff["NGC5846"]]	#Within 4arcsec
#				,["Tang, B. T., et al. 2009", 2.5/Reff['NGC5846'], 2.19, [0.,0.], 	Reff['NGC5846']]	#Within 2.5arcsec
				,['Conroy et al. (2012)', 1./8., 0.17, [0.,0.],		78.6]
				,["Gallagher et al. (2008)", 11.35/Reff['NGC5846'], -0.15, [0.10, 0.05],		Reff['NGC5846']]	#Within 11.35arcsec
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.43, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
				,["Idiart et al. (2007)", 1./8., 0.06, [0., 0.],	138.]
				,["Annibali et al. (2007)", 1./8., 0.03, [0.01, 0.01],	62.7]
			],
	"NGC5866":	[	[r"Sil'chenko (2006)", 7./Reff['NGC5866'], 0.2, [0.,0.], Reff['NGC5866']]
#				,["Terlevich, A. I., et al. 2002", 1./8., 0.44, [0., 0.],	R]	#MISSING RADIUS (IN GONZALES THESIS)
			],
	"NGC7457":	[	[r"Sil'chenko (2006)", 4./Reff['NGC7457'], 0., [0.,0.],		7./Reff['NGC7457'], -0.1, [0.,0.], Reff['NGC7457']]
				,["Serra et al. (2008)", 1./16., 0.08, [0.04,0.04],		1./2., -0.26, [0.05,0.02],	30.]
			]
}


galLiterature_metallicity_profiles ={	"NGC821": ['Proctor, R., et al. 2005', [-1.7,0.],[0.4,-0.5], 50.]	#Paper, [x1,x2],[y1,y2], 	Reff, 	Paper2, ..
					,"NGC1400": ['Spolaor, M., et al. 2008b', [-1.45, 0.1],[0.4,-0.3], 27.]
					,"NGC1407": ['Spolaor, M., et al. 2008b', [-1.8, -0.25],[0.4,-0.1], 72.]
					,"NGC3115": ['Norris et al. 2006 - Major', [-2.2, -0.05],[0.4,0.12],	93., 	'Norris et al. 2006 - Minor', [-1.8, 0.],[0.4,-0.35],	35.]
}


#Metallicity offsets for SAURON galaxies

MetOffsetSAURON = {'NGC821':[0.58939424, 0.0392181128908], 'NGC1023':[0.73287189, 0.0340970413472], 
                   'NGC2768':[0.40243949, 0.0282517407415], 'NGC2974':[0.91606397, 0.0685431997079], 
                   'NGC3377':[0.25151731, 0.0707427419238], 'NGC4278':[0.94182406, 0.0342921498452], 
                   'NGC4374':[0.92576376, 0.128961768147], 'NGC4473':[0.49067337, 0.0502709820957], 
                   'NGC4526':[0.72095879, 0.0326545489631], 'NGC5846':[1.0343796, 0.0260918621436], 
                   'NGC7457':[0.3227106, 0.0847837441226]}


#Points excluded from the dataset (indices)
dicExcluded = {'NGC1023': 27, 'NGC2974': 11, 'NGC3377': 9, 
               'NGC3607': 12, 'NGC4494': 26, 'NGC4526': 25, 'NGC4697': 42}

#Parameters kriging mapping
Theta_Kriging = {"NGC720": 10., "NGC821": 10., "NGC1023": 10.,"NGC1400":10.,
        "NGC1407":10., "NGC2768": 10., "NGC2974": 10., "NGC3115": 10.,
        "NGC3377": 10., "NGC3607":10.,"NGC3608":10., "NGC4111":10., "NGC4278": 10., 
        "NGC4365":10., "NGC4374":10., "NGC4473":10., "NGC4449":10, "NGC4486":10,  
        "NGC4494":10., "NGC4526":10., "NGC4564":10, "NGC4649":10., "NGC4697":10., 
        "NGC5846":15., "NGC5866":10., "NGC5907":10., "NGC7457":25.}


coeff_Kriging = {"NGC720": 3., "NGC821": 3., "NGC1023": 3.,"NGC1400":3.,
        "NGC1407":3., "NGC2768": 3., "NGC2974": 3., "NGC3115": 3.,
        "NGC3377": 3., "NGC3607":3., "NGC4111":3., "NGC4278": 3., 
        "NGC4365":3., "NGC4374":3., "NGC4473":3.,  
        "NGC4494":3., "NGC4526":3., "NGC4649":3., "NGC4697":3., 
        "NGC5846":3., "NGC5866":3., "NGC7457":3.}


GCsplitColour = {'NGC821': 0.97, 'NGC1400': 0.95, 'NGC1407': 0.98, 'NGC2768': 0.95, 
		'NGC3115': 0.93, 'NGC3377': 0.93, 'NGC4278': 0.86, 'NGC4365': 0.91, 
		'NGC4494': 0.99, 'NGC4649': 0.99, 'NGC5846': 0.96,	#From Usher+12
		'NGC4473': 0.99}	#Extra from Chris (source unknown)



#Cut limits in log10(R/Reff)
xLimSauron = {'NGC5846':[numpy.nan,-0.31]}




