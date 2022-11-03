"""
Module:   kPolyFitsGeom algorithm
Author:   Karunakar
Ticket:   github.com/karunakar2

Matrix of geometries from one table fitting into another table set
"""

#standard functions
import math
import random

#external modules
import numpy as np
import geopandas as gpd

#targetted imports
from shapely import affinity
from shapely.geometry import Point, LineString
from dataclasses import dataclass

#"""
class kPFGException(Exception):
    print('see kPolyFitGeom.notes() for information')
    print('no unit tests attached')
    pass
#"""

@dataclass
class kPolyFitsGeom:
    host: gpd.GeoDataFrame() = gpd.GeoDataFrame()
    tenant: gpd.GeoDataFrame() = gpd.GeoDataFrame()
    clearance: float = 0
    suppressText:bool = False

    def credits(self):
        print('Check polygons fit in another algorithm')
        print('Developed by karunakar for WDC LCM programme of work')
        print('https://github.com/karunakar2')

    def notes():
        print(15*'-')
        print('Notes:')
        print('the sub system expects geopandas geodataframes as input')

    def __del__(self):
        self.host = None
        self.tenant = None
        if self.suppressText is False:
            self.credits()

    def __enter__(self):
        #print('This module doesnot support WITH statement as of now')
        raise kPFGException('This module doesnot support WITH statement as of now')
    def __exit__(self):
        self.__del__()

    def _file_checks(self):
        try:
            if self.host.geometry.is_valid.any() == False:
                raise kPFGException('host has complex geometry')
            if self.tenant.geometry.is_valid.any() == False:
                raise kPFGException('tenant has complex geometry')
        except Exception as er:
            print(er)

    def _env_checks(self):
        #alt os.getenv('JUPYTERHUB_API_TOKEN') is None
        try:
            __IPYTHON__
            self._in_ipython_session = True
        except NameError:
            self._in_ipython_session = False
        except Exception as er:
            print(er)

    ##-------------------
        #have to consolidate functions below
    #
    #https://stackoverflow.com/a/66118219/20171247
    def _azimuth_line(self,point1, point2):
        """azimuth between 2 points (interval 0 - 180)"""
        angle = np.arctan2(point2[0] - point1[0], point2[1] - point1[1])
        return np.degrees(angle) if angle > 0 else np.degrees(angle) + 180

    """
    def _dist_line(self,point1, point2):
        #distance between points
        return math.hypot(point2[0] - point1[0], point2[1] - point1[1])
    """
    
    def _azimuth_box(self,thisRect):
        """azimuth of minimum_rotated_rectangle, angle of long side"""
        bbox = list(thisRect.exterior.coords)
        #axis1 = self._dist_line(bbox[0], bbox[3])
        #axis2 = self._dist_line(bbox[0], bbox[1])
        axis1 = Point(bbox[0]).distance(Point(bbox[3]))
        axis2 = Point(bbox[0]).distance(Point(bbox[1]))
        if axis1 <= axis2:
            az = self._azimuth_line(bbox[0], bbox[1])
        else:
            az = self._azimuth_line(bbox[0], bbox[3])
        return az

    #--- core non exposed methods
    def _preproc_gpdf(self,thisGdf,thisInnerClearance):
        thisGdf['azimuth'] = thisGdf.geometry.apply(lambda x: self._azimuth_box(x.minimum_rotated_rectangle))
        thisGdf['noOrientation'] = thisGdf.apply(lambda x: (affinity.rotate(x.geometry,-1*x.azimuth, origin='centroid')).buffer(-1*thisInnerClearance), axis=1)
        thisGdf['centerAt'] = thisGdf['noOrientation'].centroid
        return thisGdf

    def _fit_geom_check(self,thisPlot):
        tempDict=thisPlot.to_dict()
        localCopy = self.tenant.copy()
        localCopy['shiftX'] = localCopy.apply(lambda a: thisPlot.centerAt.x-a.centerAt.x, axis=1)
        localCopy['shiftY'] = localCopy.apply(lambda a: thisPlot.centerAt.y-a.centerAt.y, axis=1)
        localCopy['shiftedGeom'] = localCopy.apply(lambda x: affinity.translate(x['noOrientation'],x['shiftX'],x['shiftY']),axis=1) #or should we do plinthPlan itself: no this is simple
        localCopy['fitsHost'] = localCopy.apply(lambda x: tempDict['noOrientation'].contains(x['shiftedGeom']), axis=1)
        temp = localCopy[localCopy['fitsHost']==True]['Id'].values #.astype('int')
        if len(temp) > 0:
            return temp
        else :
            return np.nan #in case nothing fits there
    
    #--- core function exposed to run the functionality
    def _fitMatrix(self):
        return self.host[['Id','geomList']]
    matrix = property(_fitMatrix,)

    def _randomiseId(self):
        self.host['rand'] = self.host['geomList'].apply(lambda x: random.choice(list(map(int, x))))
        return self.host[['Id','rand']]
    randomId = property(_randomiseId,)

    def _randomGeom(self):
        if 'rand' not in self.host.columns:
            self._randomiseId()
        fitSet = self.host[['Id', 'azimuth', 'centerAt', 'rand']]
        fitSet = fitSet.merge(self.tenant[['Id','noOrientation', 'azimuth', 'centerAt']]
                                  ,how='inner',left_on='rand',right_on='Id', suffixes=('_h', '_t')) #h-host,t-tenant
        fitSet['shiftX'] = fitSet.apply(lambda a: a['centerAt_h'].x - a['centerAt_t'].x,axis=1)
        fitSet['shiftY'] = fitSet.apply(lambda a: a['centerAt_h'].y - a['centerAt_t'].y,axis=1)
        fitSet['shiftedGeom'] = fitSet.apply(lambda x: affinity.translate(x['noOrientation'],x['shiftX'],x['shiftY']),axis=1)
        fitSet['alignedGeom'] = fitSet.apply(lambda x: affinity.rotate(x.shiftedGeom, x.azimuth_h, origin='centroid'), axis=1)
        tempDf = fitSet[['Id_h','Id_t','alignedGeom']].copy()
        return gpd.GeoDataFrame(tempDf,geometry='alignedGeom')
    sampleFittedGeometry = property(_randomGeom,)
                          
    def run(self):
        self._env_checks()
        self._file_checks()

        #display progress module
        if self._in_ipython_session: #jupyter notebook
            from tqdm.autonotebook import tqdm
        else:
            from tqdm.auto import tqdm
        tqdm.pandas()

        self._preproc_gpdf(self.host,self.clearance)
        self.tenant = self.tenant[self.tenant.geometry.area <= self.host.geometry.area.max()]
        self._preproc_gpdf(self.tenant,0)
        
        #print(self.host.head(),self.tenant.head())
        self.host['geomList'] = self.host.progress_apply(lambda x: self._fit_geom_check(x), axis=1)
        self.host= self.host.dropna()

        #print(self.host.head())
        print('finished, post process...')

