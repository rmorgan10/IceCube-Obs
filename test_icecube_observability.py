import os

"""
# TEST
event_dict = {'80305071': {'ra': 98.3268, # EHE
                           'dec': -14.4861, 
                           'eventid': 80305071, 
                           'time': None},
              '80127519': {'ra': 45.8549,
                           'dec': 15.7851, 
                           'eventid': 80127519, 
                           'time': '2016/12/10 20:07:16'},
              '26552458': {'ra': 122.7980,
                           'dec': -0.7331, 
                           'eventid': 26552458, 
                           'time': None},
              '6888376': {'ra': 214.5440,
                          'dec': -0.3347, 
                          'eventid': 6888376, 
                          'time': None},
              '32674593': {'ra': 221.6750, # HESE
                           'dec': -26.0359, 
                           'eventid': 32674593, 
                           'time': None},
              '65274589': {'ra': 304.7300,
                           'dec': -26.2380, 
                           'eventid': 65274589, 
                           'time': None},
              '38561326': {'ra': 40.8252,
                           'dec': 12.5592, 
                           'eventid': 38561326, 
                           'time': None},
              '58537957': {'ra': 199.3100,
                           'dec': -32.0165, 
                           'eventid': 58537957, 
                           'time': None},
              #'6888376': {'ra': 215.1090,
              #            'dec': -0.4581, 
              #            'eventid': 6888376, 
              #            'time': None},
              '67093193': {'ra': 240.5683,
                           'dec': 9.3417, 
                           'eventid': 67093193, 
                           'time': None}}
"""
"""
              '999': {'ra': 45.8549 + 180.,
                      'dec': 15.7851, 
                      'eventid': 999, 
                      'time': '2016/12/10 20:07:16'}}
"""
# ACTUAL TRIGGER
#event_dict = {'50579430': {'ra': 77.2853,
#                           'dec': 5.7517, 
#                           'eventid': 50579430, 
#                           'time': None}}
# Updated: using best-fit location 
#event_dict = {'50579430': {'ra': 77.43,
#                           'dec': 5.72, 
#                           'eventid': 505794,
#                           'time': None}}
# Updated: trying to fully contain the 90% ellipse within DECam field of view
#event_dict = {'50579430': {'ra': 77.68,
#                           'dec': 5.87, 
#                           'eventid': 505794,
#                           'time': None}}
#event_dict = {'56068624': {'ra': 162.5790,
#                           'dec': -15.8611, 
#                           'eventid': 56068624, 
#                           'time': None}}
# Initial coordinates
#event_dict = {'17569642': {'ra': 340.2500,
#                           'dec': 7.3140, 
#                           'eventid': 17569642, 
#                           'time': None}}
# Updated coordinates from offline reconstruction
#event_dict = {'17569642': {'ra': 340.00,
#                           'dec': 7.40, 
#                           'eventid': 17569642, 
#                           'time': None}}
#event_dict = {'17569642': {'ra': 340.00,
#                           'dec': 7.40, 
#                           'eventid': 17569642, 
#                           'time': '2017/11/29 00:00:00'}}
# Event 2018 Sep 8
#event_dict = {'15947448': {'ra': 337.68,
#                           'dec': -20.70, 
#                           'eventid': 15947448, 
#                           'time': '2019/3/31 00:27:07'}}

#GW
#event_dict = {'190510': {'ra': 85.0,
#                           'dec': -30.,
#                           'eventid': 190510,
#                           'time': '2019/5/11 00:27:07'}}

#event_dict = {'42419327': {'ra': 120.304,
#                           'dec': 6.3568,
#                           'eventid': 42419327,
#                           'time': None}}

event_dict = {'42419327': {'ra': 342.7772,
                           'dec': 10.0547,
                           'eventid': 132707,
                           'time': None}}

#IC190504
#event_dict = {'766165': {'ra': 65.7866,
#                         'dec': -37.4431,
#                         'eventid': 766165,
#                         'time': '2019/5/25 00:27:07'}}

#print event_dict.keys()
print(event_dict.keys())

for key in event_dict.keys():
    #if key != '999':
    #    continue

    if event_dict[key]['time'] is not None:
        command = 'pythonw icecube_observability.py --ra %.4f --dec %.4f --eventid %i --time \"%s\"'%(event_dict[key]['ra'],
                                                                                                        event_dict[key]['dec'],
                                                                                                        event_dict[key]['eventid'],
                                                                                                        event_dict[key]['time'])
    else:
        #print "Let's begin"
        print("Let's begin")
        command = 'python icecube_observability.py --ra %.4f --dec %.4f --eventid %i'%(event_dict[key]['ra'],
                                                                                          event_dict[key]['dec'],
                                                                                          event_dict[key]['eventid'])

    os.system(command)
