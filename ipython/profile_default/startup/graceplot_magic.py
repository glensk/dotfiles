import numpy as np
import GracePlot as gp

#x=np.array([1,2,3,4,5,6,7,8,9,10])
#y=np.array([9,2,1,1,5,8,7,8,9,10])


def plotxy(x,y):
    #xplot = gp.GracePlot() # A grace session opens
    gp.GracePlot().plot(gp.Data(x=x,y=y))

def ploty(y):
    #xplot = gp.GracePlot() # A grace session opens
    gp.GracePlot().plot(gp.Data(x=np.arange(len(y)),y=y))

#plotx(x,y)

