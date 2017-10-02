
# coding: utf-8

# In[233]:


import sys
from random import choice

#Stations
station = 'SYD|-33.86,151.21,39 MEL|-37.81,144.96,119 ADE|-34.92,138.60,8 BRB|-27.23,153.8,10 DAR|-12.25,130.54,35 QUN|-17.39,141.5,9 CAN|-35.19,149.12,578 WAU|-35.2,117.55,13 PER|-31.56,115.59,20 CAL|-16.52,145.45,8'.split()
longitude = [-33.86,-37.81,-34.92,-27.23,-12.25,-17.39,-35.19,-35.2,-31.56,-16.52]

#Local Time
day = range(1, 32)
month = range(1, 13)
spring_autumen, summer, winter = [3, 4, 9, 10], [5, 6, 7, 8], [11, 12, 1 ,2]
year = [2016, 2017]

seconds = range(0, 60)
minute = range(0, 60)
hour = range(24)

#Conditions
conditions = 'Rain Snow Sunny'.split()

#Temperature
spring_autumn_temp, summer_temp, winter_temp = range(17, 30), range(20, 50), range(-3, 17)

#Pressure
pressure = range(100, 1200)

#Humidity
humidity = range(30, 97)
   

    
import math
from math import pi

import calendar

class Sun:
    
    def __init__(self):

	
	# Some conversion factors between radians and degrees
	self.RADEG = 180.0 / pi
	self.DEGRAD = pi / 180.0
	self.INV360 = 1.0 / 360.0
	

    def daysSince2000Jan0(self, y, m, d):

	return (367*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530)
    

    # The trigonometric functions in degrees
    def sind(self, x):

	return math.sin(x * self.DEGRAD)

    def cosd(self, x):

	return math.cos(x * self.DEGRAD)

    def tand(self, x):

	return math.tan(x * self.DEGRAD)

    def atand(self, x):

	return math.atan(x) * self.RADEG
    
    def asind(self, x):

	return math.asin(x) * self.RADEG

    def acosd(self, x):

	return math.acos(x) * self.RADEG

    def atan2d(self, y, x):

	return math.atan2(y, x) * self.RADEG

                            

    def dayLength(self, year, month, day, lon, lat):

	return self.__daylen__(year, month, day, lon, lat, -35.0/60.0, 1)
    
    def dayCivilTwilightLength(self, year, month, day, lon, lat):

	return self.__daylen__(year, month, day, lon, lat, -6.0, 0)

    def dayNauticalTwilightLength(self, year, month, day, lon, lat):

	return self.__daylen__(year, month, day, lon, lat, -12.0, 0)

    def dayAstronomicalTwilightLength(self, year, month, day, lon, lat):

	return self.__daylen__(year, month, day, lon, lat, -18.0, 0)
    
    def sunRiseSet(self, year, month, day, lon, lat):

	return self.__sunriset__(year, month, day, lon, lat, -35.0/60.0, 1)

    
    def civilTwilight(self, year, month, day, lon, lat):

	return self.__sunriset__(year, month, day, lon, lat, -6.0, 0)
    
    def nauticalTwilight(self, year, month, day, lon, lat):

	return self.__sunriset__(year, month, day, lon, lat, -12.0, 0)
    
    def astronomicalTwilight(self, year, month, day, lon, lat):

	return self.__sunriset__(year, month, day, lon, lat, -18.0, 0)
    
    # The "workhorse" function for sun rise/set times
    def __sunriset__(self, year, month, day, lon, lat, altit, upper_limb):

	# Compute d of 12h local mean solar time
	d = self.daysSince2000Jan0(year,month,day) + 0.5 - (lon/360.0)
	
	# Compute local sidereal time of this moment 
	sidtime = self.revolution(self.GMST0(d) + 180.0 + lon)
	
	# Compute Sun's RA + Decl at this moment 
	res = self.sunRADec(d)
	sRA = res[0]
	sdec = res[1]
	sr = res[2]
	
	# Compute time when Sun is at south - in hours UT 
	tsouth = 12.0 - self.rev180(sidtime - sRA)/15.0;
	
	# Compute the Sun's apparent radius, degrees 
	sradius = 0.2666 / sr;
	
	# Do correction to upper limb, if necessary 
	if upper_limb:
	    altit = altit - sradius

	
	cost = (self.sind(altit) - self.sind(lat) * self.sind(sdec))/               (self.cosd(lat) * self.cosd(sdec))

	if cost >= 1.0:
	    rc = -1 
	    t = 0.0           # Sun always below altit
	    
	elif cost <= -1.0:
	    rc = +1 
	    t = 12.0;         # Sun always above altit

	else:
	    t = self.acosd(cost)/15.0   # The diurnal arc, hours

	
	# Store rise and set times - in hours UT 
	return (tsouth-t, tsouth+t)


    def __daylen__(self, year, month, day, lon, lat, altit, upper_limb):

	
	# Compute d of 12h local mean solar time 
	d = self.daysSince2000Jan0(year,month,day) + 0.5 - (lon/360.0)
		 
	# Compute obliquity of ecliptic (inclination of Earth's axis) 
	obl_ecl = 23.4393 - 3.563E-7 * d
	
	# Compute Sun's position 
	res = self.sunpos(d)
	slon = res[0]
	sr = res[1]
	
	# Compute sine and cosine of Sun's declination 
	sin_sdecl = self.sind(obl_ecl) * self.sind(slon)
	cos_sdecl = math.sqrt(1.0 - sin_sdecl * sin_sdecl)
	
	# Compute the Sun's apparent radius, degrees 
	sradius = 0.2666 / sr
	
	# Do correction to upper limb, if necessary 
	if upper_limb:
	    altit = altit - sradius

	    
	cost = (self.sind(altit) - self.sind(lat) * sin_sdecl) /                (self.cosd(lat) * cos_sdecl)
	if cost >= 1.0:
	    t = 0.0             # Sun always below altit
	
	elif cost <= -1.0:
	    t = 24.0      # Sun always above altit
	
	else:
	    t = (2.0/15.0) * self.acosd(cost);     # The diurnal arc, hours
	    
	return t

    
    def sunpos(self, d):
 
	M = self.revolution(356.0470 + 0.9856002585 * d)
	w = 282.9404 + 4.70935E-5 * d
	e = 0.016709 - 1.151E-9 * d
	
	# Compute true longitude and radius vector 
	E = M + e * self.RADEG * self.sind(M) * (1.0 + e * self.cosd(M))
	x = self.cosd(E) - e
	y = math.sqrt(1.0 - e*e) * self.sind(E)
	r = math.sqrt(x*x + y*y)              #Solar distance 
	v = self.atan2d(y, x)                 # True anomaly 
	lon = v + w                        # True solar longitude 
	if lon >= 360.0:
	    lon = lon - 360.0   # Make it 0..360 degrees
	    
	return (lon,r)
    

    def sunRADec(self, d):
 
	res = self.sunpos(d)
	lon = res[0]  # True solar longitude
	r = res[1]    # Solar distance
	
	# Compute ecliptic rectangular coordinates (z=0) 
	x = r * self.cosd(lon)
	y = r * self.sind(lon)
	
	# Compute obliquity of ecliptic (inclination of Earth's axis) 
	obl_ecl = 23.4393 - 3.563E-7 * d
	
	# Convert to equatorial rectangular coordinates - x is unchanged 
	z = y * self.sind(obl_ecl)
	y = y * self.cosd(obl_ecl)

	# Convert to spherical coordinates 
	RA = self.atan2d(y, x)
	dec = self.atan2d(z, math.sqrt(x*x + y*y))

	return (RA, dec, r)


    def revolution(self, x):

	return (x - 360.0 * math.floor(x * self.INV360))
    
    def rev180(self, x):
	
	return (x - 360.0 * math.floor(x * self.INV360 + 0.5))

    def GMST0(self, d):

	# Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  
	# L = M + w, as defined in sunpos().  Since I'm too lazy to 
	# add these numbers, I'll let the C compiler do it for me.  
	# Any decent C compiler will add the constants at compile   
	# time, imposing no runtime or code overhead.               
						
	sidtim0 = self.revolution((180.0 + 356.0470 + 282.9404) +
	                             (0.9856002585 + 4.70935E-5) * d)
	return sidtim0;

    def solar_altitude(self, latitude, year, month, day):
   
        N = self.daysSince2000Jan0(year, month, day)
        res =  self.sunRADec(N)
        declination = res[1]
        sr = res[2]

        # Compute the altitude
        altitude = 90.0 - latitude  + declination

        # In the tropical and  in extreme latitude, values over 90 may occurs.
        if altitude > 90:
            altitude = 90 - (altitude-90)

        if altitude < 0:
            altitude = 0

        return altitude

    def get_max_solar_flux(self, latitude, year, month, day):
     
        (fEot, fR0r, tDeclsc) = self.equation_of_time(year, month, day, latitude)
        fSF = (tDeclsc[0]+tDeclsc[1])*fR0r

        # In the case of a negative declinaison, solar flux is null
        if fSF < 0:
            fCoeff = 0
        else:
            fCoeff =  -1.56e-12*fSF**4 + 5.972e-9*fSF**3 -                     8.364e-6*fSF**2  + 5.183e-3*fSF - 0.435
       
        fSFT = fSF * fCoeff 

        if fSFT < 0:
            fSFT=0

        return fSFT

    def equation_of_time(self, year, month, day, latitude):
      
        nJulianDate = self.Julian(year, month, day)
        # Check if it is a leap year
        if(calendar.isleap(year)):
            fDivide = 366.0
        else:
            fDivide = 365.0
        # Correction for "equation of time"
        fA = nJulianDate/fDivide*2*pi
        fR0r = self.__Solcons(fA)*0.1367e4
        fRdecl = 0.412*math.cos((nJulianDate+10.0)*2.0*pi/fDivide-pi)
        fDeclsc1 = self.sind(latitude)*math.sin(fRdecl)
        fDeclsc2 = self.cosd(latitude)*math.cos(fRdecl)
        tDeclsc = (fDeclsc1, fDeclsc2)
        # in minutes
        fEot = 0.002733 -7.343*math.sin(fA)+ .5519*math.cos(fA)                - 9.47*math.sin(2.0*fA) - 3.02*math.cos(2.0*fA)                - 0.3289*math.sin(3.*fA) -0.07581*math.cos(3.0*fA)                -0.1935*math.sin(4.0*fA) -0.1245*math.cos(4.0*fA)
        # Express in fraction of hour
        fEot = fEot/60.0
        # Express in radians
        fEot = fEot*15*pi/180.0

        return (fEot, fR0r, tDeclsc)

    def __Solcons(self, dAlf):
       
        
        dVar = 1.0/(1.0-9.464e-4*math.sin(dAlf)-0.01671*math.cos(dAlf)-                     + 1.489e-4*math.cos(2.0*dAlf)-2.917e-5*math.sin(3.0*dAlf)-                     + 3.438e-4*math.cos(4.0*dAlf))**2
        return dVar


    def Julian(self, year, month, day):
       
        if calendar.isleap(year): # Bissextil year, 366 days
            lMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        else: # Normal year, 365 days
            lMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

        nJulian = lMonth[month-1] + day
        return nJulian



if __name__ == "__main__":

    k = Sun()
    print k.get_max_solar_flux(-34.08, 2017, 10, 01)
    
          
    
    


# In[236]:

import random
def stationmapping(x):
    return {
        -33.86: 'SYD|-33.86,151.21,39',
        -37.81: 'MEL|-37.81,144.96,119',
        -34.92: 'ADE|-34.92,138.60,8',
        -27.23: 'BRB|-27.23,153.8,10',
        -12.25: 'DAR|-12.25,130.54,35',
        -17.39: 'QUN|-17.39,141.5,9',
        -35.19: 'CAN|-35.19,149.12,578',
        -35.2 : 'WAU|-35.2,117.55,13',
        -31.56: 'PER|-31.56,115.59,20',
        -16.52: 'CAL|-16.52,145.45,8',
    }[x]

def generate_data(i):
    mnt = choice(month)
    if mnt == 2:
        dy = choice(day)
        while dy > 28:
            dy = choice(day)
    else:
        dy = choice(day)
    yr = choice(year)

    sec = choice(seconds)
    min = choice(minute)
    hr = choice(hour)

    cnd = choice(conditions)
    hmd = choice(humidity)
    temp = 30
    #local_time = '%s-%s-%sT%s:%s:%sZ' %(yr, mnt.zfill(2), dy.zfill(2), hr.zfill(2), min.zfill(2), sec.zfill(2))
    solar2 = []

#    for i in range(10):
    solar3 = (k.solar_altitude(longitude[i], yr, mnt, dy))  
    #solar2.append(solar3) #Set temperature as per seasons
    if int(solar3) in range(47,67):
        temp = choice(spring_autumn_temp)
            
    elif int(solar3) in range(67,90):
        #while cnd == 'Snow':
        #    cnd = choice(conditions)
        temp = choice(summer_temp)
    elif int(solar3) in range(20,47):
        temp = choice(winter_temp)
    #Drop temp based on weather conditions
    if cnd == 'Snow':
        #temp = 0
        temp -= choice(range(2, 7))
    elif cnd == 'Rain':
        temp -= choice([1, 2, 3])
        hmd += choice(range(10, 20))
        #Set humidity to max(97) if it crosses 97
        hmd = max(97, hmd)
    # Drop night temeprature by 3 to 7 degrees
    if hr in range(8, 22):
        temp += random.choice(range(2, 5))
    dy, mnt, yr, sec, min, hr = map(str, [dy, mnt, yr, sec, min, hr])
    local_time = '%s-%s-%sT%s:%s:%sZ' %(yr, mnt.zfill(2), dy.zfill(2), hr.zfill(2), min.zfill(2), sec.zfill(2))  
    st1 = stationmapping(longitude[i])
    return '|'.join(map(str, [st1, local_time, cnd, temp, choice(pressure), hmd]))
       # print ('|'.join(map(str, [st1, local_time, cnd, temp, choice(pressure), hmd, solar2[i]])))

       


for i in range(int(10)):
      print generate_data(i)


# In[ ]:



