ellipsoid_a=6378137.
ellipsoid_b=6356752.
R=6371000.8
dlon=0.01   # decimal degrees
lat=37.9804   # decimal degrees
dx(lat,dlon)=R*dlon*pi/180.*cos(lat*pi/180.)
dxEq(lat,dlon)=ellipsoid_a*dlon*pi/180.*cos(lat*pi/180.)
print "lat_deg=",lat," dlon_deg=",dlon," dx=",dx(lat,dlon)," dxEq=",dxEq(lat,dlon)


