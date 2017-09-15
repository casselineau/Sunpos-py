import numpy as N

def sunpos(iYear, iMonth, iDay, dHours, dMinutes, dSeconds, dLatitude, dLongitude):
	'''
	From the C++ version at http://www.psa.es/sdg/sunpos.htm
	Paper: https://doi.org/10.1016/S0038-092X(00)00156-0
	Azimuth is clockwise from North
	'''	
	# Constants:
	dEarthMeanRadius = 6371.01	# In km
	dAstronomicalUnit = 149597890
	rad = N.pi/180.

	# Calculate difference in days between the current Julian Day 
	# and JD 2451545.0, which is noon 1 January 2000 Universal Time
	# Calculate time of the day in UT decimal hours
	dDecimalHours = dHours + (dMinutes + dSeconds / 60.0 ) / 60.0
	# Calculate current Julian Day
	liAux1 =(iMonth-14)/12
	liAux2=(1461*(iYear + 4800 + liAux1))/4 + (367*(iMonth - 2-12*liAux1))/12- (3*((iYear + 4900 + liAux1)/100))/4+iDay-32075
	dJulianDate=float(liAux2)-0.5+dDecimalHours/24.0

	# Calculate difference between current Julian Day and JD 2451545.0 
	dElapsedJulianDays = dJulianDate-2451545.0

	# Calculate ecliptic coordinates (ecliptic longitude and obliquity of the 
	# ecliptic in radians but without limiting the angle to be less than 2*Pi 
	# (i.e., the result may be greater than 2*Pi)
	dOmega=2.1429-0.0010394594*dElapsedJulianDays
	dMeanLongitude = 4.8950630+ 0.017202791698*dElapsedJulianDays # Radians
	dMeanAnomaly = 6.2400600+ 0.0172019699*dElapsedJulianDays
	dEclipticLongitude = dMeanLongitude + 0.03341607*N.sin(dMeanAnomaly) + 0.00034894*N.sin(2*dMeanAnomaly)-0.0001134 - 0.0000203*N.sin(dOmega)
	dEclipticObliquity = 0.4090928 - 6.2140e-9*dElapsedJulianDays + 0.0000396*N.cos(dOmega)

	# Calculate celestial coordinates ( right ascension and declination ) in radians 
	# but without limiting the angle to be less than 2*Pi (i.e., the result may be 
	# greater than 2*Pi)
	dSin_EclipticLongitude = N.sin(dEclipticLongitude)
	dY = N.cos(dEclipticObliquity) * dSin_EclipticLongitude
	dX = N.cos(dEclipticLongitude)
	dRightAscension = N.arctan2(dY, dX)
	if (dRightAscension<0.0):
		dRightAscension = dRightAscension + 2*N.pi
	dDeclination = N.arcsin(N.sin(dEclipticObliquity)*dSin_EclipticLongitude)

	# Calculate local coordinates ( azimuth and zenith angle ) in degrees
	dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283*dElapsedJulianDays + dDecimalHours
	dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime*15 + dLongitude)*rad
	dHourAngle = dLocalMeanSiderealTime - dRightAscension
	dLatitudeInRadians = dLatitude*rad
	dCos_Latitude = N.cos( dLatitudeInRadians )
	dSin_Latitude = N.sin( dLatitudeInRadians )
	dCos_HourAngle= N.cos( dHourAngle )
	dZenithAngle = (N.arccos(dCos_Latitude*dCos_HourAngle*N.cos(dDeclination) + N.sin(dDeclination)*dSin_Latitude))
	dY = -N.sin( dHourAngle )
	dX = N.tan(dDeclination)*dCos_Latitude - dSin_Latitude*dCos_HourAngle
	dAzimuth = N.arctan2(dY, dX)
	if (dAzimuth<0.0):
		dAzimuth = dAzimuth + 2*N.pi

	# Parallax Correction
	dParallax=(dEarthMeanRadius/dAstronomicalUnit)*N.sin(dZenithAngle)
	dZenithAngle = dZenithAngle + dParallax

	return dAzimuth, dZenithAngle # in rads

