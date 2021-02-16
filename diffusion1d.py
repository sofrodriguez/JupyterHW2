def oso(D,age):
# outputs: 
# topo -- topographic elevation of the cone [m]
# distance -- distance along the cone profile [m]
# coneslope -- mean slope of cone [deg]
# inputs: 
# D -- landscape diffusivity [m^2/kyr]
# age -- age of the landform (time of simulation) [kyr]

    import numpy as np

    # Topo and Grid Parameters
    dx=3
    topo=np.loadtxt('osotransect.txt')
    #opo=np.loadtxt('rowantransect.txt')
    nx=len(topo)
    
    
    distance=np.zeros(nx)
    distance=np.arange(0,nx*dx,dx).tolist()

    # 1d Diffusion
    tend=age           # kyr

    toponew=np.zeros(nx)
    dt=0.1*dx*dx/D
    slope=np.zeros(nx)

    t=0
    while t<tend:
    
        for i in range(1,nx-1):
            toponew[i]=topo[i]+D*dt/(dx*dx)*(topo[i+1]-2*topo[i]+topo[i-1])

        toponew[0]=topo[0]+D*dt/(dx*dx)*(topo[1]-2*topo[0]+topo[0])
        toponew[nx-1]=topo[nx-1]+D*dt/(dx*dx)*(topo[nx-1]-2*topo[nx-1]+topo[nx-2])
    
        t+=dt
        topo=toponew[:]
        
    slope=np.gradient(topo,dx)
    slope=180/np.pi*np.arctan(slope)


    # Compute standard deviation of slope over 30m interval
    interval=10
    SDS=np.array([])
    for i in range(interval,len(topo),interval):
        SDS=np.concatenate((SDS,[np.std(slope[i-interval:i], ddof=1)]))
        
    SDS=np.mean(SDS)        

    return (distance,topo,SDS);

