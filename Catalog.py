class Catalog():
    
    def __init__(self, df, num, time = np.linspace(0, 2*np.pi, 768), freqs=np.linspace(150, 160, 10)):
        """
        Initializes the Catalog and makes a new Pandas DataFrame containing 
        the Right Ascension and Declination as well as all the fluxes
        
        Parameters:
        df = Pandas DataFrame
        num = sample size num
        time = time range array
        freqs= freq array
        Creates: 
        self.orgFrame as Pandas DataFrame type, cleaned, all orginal data, with size num
        self.flux as np array type, cleaned, contains only sample size num and freq freqs
        self.RA = right ascencion from J2000 in degrees as numpy.ndarray
        self.DEC = declination from J2000 in degrees as numpy.ndarray
        self.orgFreqMHz = all old freq in order as in df passed in
        self.FreqMHz = array of all freq run
        self.time = array of all time run
        self.altaz = [alt, az] for each source for each time
        self.kAlphas = first row is k and second row is alpha for all sources
        
        Functions:
        get_flux(self, freq) return an array of flux at Freq
        getAltazAtFreq(self, freq) return an array of [alt, az] at Freq
        kAlpha(self) return a size 2 array with first being k of all sources and second being alphas
        position_vectors()
        """
        d = df.sort_values(by = "int_flux_076", ascending=False)
        d = d.dropna(axis = 0, how = "any")
        self.orgFrame = d[[col for col in d.keys() if 
                           ((col=='RAJ2000' or col=='DEJ2000'or col[:9]=='int_flux_') 
                            and col!= 'int_flux_fit_200') and col != 'int_flux_wide']][:num]
        self.orgFreqMHz = [(int) (col[10:]) for col in self.orgFrame.keys() if (col[:9]=='int_flux_')]
        self.RA = np.asarray(self.orgFrame.RAJ2000)
        self.DEC = np.asarray(self.orgFrame.DEJ2000)
        self.time = time
        print(time.shape)
        self.alt_az = self.altaz(time)
        self.k, self.alpha = self.kAlpha()
        self.freqMHz = freqs
        self.s_vectors = self.position_vector()
        
        self.flux = self.get_flux()
        
        
    def altaz(self, time, rad = False, lat = 37.875*np.pi/180):
        """
        modified from from Max
        Calculates Altitude and Azimuth at given times, centered at HERA
        parameters: 
        time = lst in radian
        ra = Right ascension of the star in degrees or radians
        dec = Declination of the star in degrees or radians
        rad = True if RA and DEC are in radian or False if RA and DEC are in degree
        returns:
        Returns alt, and az all in radian
        """
        ra = self.RA
        dec = self.DEC
        PI = np.pi
        alt = []
        az = []
        #Setting the latitude/longitude of HERA
        if rad == False:
            ra =ra * PI/180
            dec =dec * PI/180
        for t in time:
            #Converts the LST to hour angle
            hour_array = np.mod((t - ra), 2*PI)

            #Calculates the Altitude and Azimuth
            alt.append( np.arcsin(np.sin(dec)*np.sin(lat)
                            + np.cos(dec)*np.cos(lat)*np.cos(hour_array)))
            az.append(np.arctan2(np.sin(hour_array)*np.cos(dec),
                            np.cos(hour_array)*np.cos(dec)*np.sin(lat) - np.sin(dec)*np.cos(lat)) + PI)
        
        array = []
        self.alt=np.array(alt)
        self.az=np.array(az)
        for alti, azi in zip(self.alt, self.az):
            array.append([alti, azi])
        return array
    
    def position_vector(self):
        """
        Calculates the position vector to a source in the sky 
        Parameters: 
        alt = the altitude of the source from the horizon, given in radians            
        az  = the angle from north, moving towards the east, given in radians
        Returns
        vector in cartesian coordinates [x,y,z] as a numpy.ndarray
        """
        x = np.cos(self.alt) * np.cos(self.az)
        y = np.cos(self.alt) * np.sin(self.az)
        z = np.sin(self.alt) 
        ##vector = np.array([x,y,z])
        return np.array([x, y, z])
    
    def kAlpha(self):
        """
        return an array of k and alphas
        first column is k for all sources and second column is a for all sources
        """
        flux1 = np.array(self.orgFrame['int_flux_076'])
        flux2 = np.array(self.orgFrame['int_flux_084'])
        k = []
        alpha = []
        a = np.log(flux1/flux2)/np.log(84e6/76e6)
        k = flux1 / np.power(76, -a)
        return k, a
    
    def get_flux(self):
        return np.einsum('i, ki->ik', self.k, np.array([np.power(freq, -self.alpha) for freq in self.freqMHz]))
    
    def getPhase(b):
        """return an array phase factors with number of row being number of source, 
        number of column being number of freq"""       
        arr_source = []
        for s in self.s_vectors:
            arr_freq = []
            for freq in self.freqMHz:
                arr_freq.append(np.exp(-2*np.pi*1.j*freq/spc.c*np.dot(s, b)))
            arr_source.append(arr_freq)
        return arr_source