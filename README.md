The package has several external module dependencies.
This module has been tested for basic functionality and reproduces reasonable results
future version to include tests

The package implementation is straight forward.
The following code is self explanatory.

Two geopandas dataframe are required.
1. the polygons table into which others are to be checked and moved.
2. the polygon that will go in and sit in the above.

from kPolyFitsGeom import kPolyFitsGeom
myAlg = kPolyFitsGeom(plotsOfInterest,homesTemplate)
myAlg.run()

def myPlots():
    ax = plotsOfInterest.plot()
    myAlg.sampleFittedGeometry.plot(ax=ax,color='k',alpha=0.33)
    plt.show()
myPlots()