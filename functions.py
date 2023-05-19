
# optimisation
def optimisation(month):
    '''
    water and energy optimisation
    input kth month to find optimised fluxes for month k
    months 1 to 144
    '''
    i = month
    # error 
    S = [Perror[i], LEerror[i], S_error[i],
         DSRerror[i], DLRerror[i], USWerror[i], ULWerror[i],
         SHerror[i],E_error[i]]
    
    S = np.square(S) 
    S_obs = np.diag(S)
    # observation vector 
    F =[Pobs[i],LEobs[i],D[i], DSRobs[i], DLRobs[i],
        USWobs[i], ULWobs[i], SHobs[i],E[i]]
    
    # budget constraint
    AT= np.array([1,-1,-1, 0, 0, 0, 0, 0, 0])
    A = np.reshape(AT,(1,9))
    
    BT= np.array([0,-1,0, 1, 1, -1, -1, -1, -1])
    B = np.reshape(BT,(1,9))

    S_Robs_inv = np.linalg.inv(S_obs)

    alpha1 = S_Robs_inv@F
    alpha2 = np.append(alpha1,0)
    alpha = np.append(alpha2,0)

    x2 = np.append(A,0)
    x2 = np.append(x2,0).reshape(11,1)
    
    x3 =np.append(B,0)
    x3 =np.append(x3,0).reshape(11,1)

    X = np.vstack((S_Robs_inv,A))
    X = np.vstack((X,B))
    
    X = np.hstack((X,x2))
    X = np.hstack((X,x3))

    X_trans = np.transpose(X)
    beta = (np.linalg.inv((X@X_trans))@X_trans)@alpha
    return beta, beta[2], beta[8]


def endo_energy_fluxesEXT(MASK):
    w2m = 28.9 # watts per square meter to mm/day
    
    DSR = mask_and_weight(DSR_global[22:,:,:], MASK)/w2m
    DLR = mask_and_weight(DLR_global[22:,:,:], MASK)/w2m
    USW = mask_and_weight(USW_global[22:,:,:], MASK)/w2m
    ULW = mask_and_weight(ULW_global[22:,:,:], MASK)/w2m
    LEm  =  mask_and_weight(LE[12:156], MASK)
    SHm =  mask_and_weight(SH[12:156], MASK)  
    
    # Mean seasonal cycle LE 
    msc_LE = monthly_mean(LEm)
     # need to add 7 more years 
    LE2 = np.array(list(msc_LE)*7)
    # add additional 7 years 
    LE2m = np.append(LEm, LE2[:-1])
    
    # Mean seasonal cycle LE 
    msc_SH = monthly_mean(SHm)
     # need to add 7 more years 
    SH2 = np.array(list(msc_SH)*7)
    # add additional 7 years 
    SH2m = np.append(SHm, SH2[:-1])



    NET = DSR[:-6] + DLR[:-6] -USW[:-6] - ULW[:-6] - LE2m[0:216] - SH2m[0:216] 
    return  DSR[:216], DLR[0:216], USW[0:216], ULW[0:216], SH2m[0:216], NET  

def endo_fluxes_EXT(MASK):

    Pm = mask_and_weight(P[12:-5], MASK) # P 2002 - 2021 227 months
    dSm =  mask_and_weight(dS[:], MASK) # dS JAN 2002 to 2020
    LEm =  mask_and_weight(LE[12:156], MASK) # dS JAN 2002 to 2020

    # Mean seasonal cycle LE 
    msc_LE = monthly_mean(LEm)
     # need to add 7 more years 
    LE2 = np.array(list(msc_LE)*7)
    # add additional 7 years 
    LE2m = np.append(LEm, LE2[:-1])
    
    Res = Pm[:] - LE2m[:]- dSm[:] #2002 to 2013
    return Pm[0:216], LE2m[0:216], dSm[0:216], Res[0:216]



  def deseason(storage):
    '''
    deseasonalised storage cycle for 18 years 
    '''
    mscSTOR = monthly_mean(storage[0:216])
    stor_deseason = storage[0:216] - np.array(list(mscSTOR)*18)[:216]
    return stor_deseason
  
  def FISd(F):
    '''
    input: flux (in units mm/day)
    out put: detrended flux inferred storage (in units cm)
    '''
    Fd = F - np.mean(F) 
    stor = np.insert(Fd[0:215]*3.046, 0, som_stor2020[0])
    FIS = np.cumsum(stor)
    return FIS

  
  def FISeYC(NET):
    '''
    input: NET energy flux (in mm/day) 
    output: flux-inferred energy storage with yearly balance
    '''
    Yearlymean = np.zeros(18)
    for i in range(18):
        Yearlymean[i] = np.mean(NET[12*i:i*12+12])
    Ym = np.repeat(Yearlymean,12)   
    NET_detrended = NET[0:216] - Ym
    X = np.insert(NET_detrended[:-1]*3.046,0,0)
    FISe = np.cumsum(X)
    return FISe

  
def som_stor2020(storage):
  '''
  input: GRACE Storage 2001-2021
  output: start of month GRACE storage 2002 - 2020
  '''
    storage2 = storage[11:240]
    som_stor = np.zeros(227)
    for i in range(227):
        som_stor[i] = ((storage2[i+1] + storage2[i])/2)
    return som_stor[0:216]

def weight_ENDO(masked_data, Region):
    '''
    weights dats according to gridbox size
    '''
    weight = mask_ENDO2D(clat2D, Region)
    a = len(masked_data[:,0,0])
    weighted_data = np.zeros(a)
    for i in range(a):
        weighted_data[i] =(masked_data[i,:,:]*weight).sum()/(weight).sum()
    return weighted_data


def mask_ENDO2D(data,Region):
    '''
    input:
    '''
    mask_boolean = np.zeros(shape =(360,720))
    for i in range(360):
        for j in range(720):
            if Region[i,j] > 0:
                mask_boolean[i,j] = 0
            else:
                mask_boolean[i,j] = 1
    masked_data = ma.masked_array(data, mask = mask_boolean)
    return masked_data

def mask_ENDO(data, Region):
    '''
    input:
    '''
    a = len(data[:,1])
    mask_boolean = np.zeros(shape =(a,360,720))
    for i in range(360):
        for j in range(720):
            if Region[i,j] > 0:
                mask_boolean[:,i,j] = 0
            else:
                mask_boolean[:,i,j] = 1
    masked_data = ma.masked_array(data, mask = mask_boolean)
    return masked_data



def mask_and_weight(data, Region):
    '''
    weights dats according to gridbox size
    '''
    a = len(data[:,1])
    mask_boolean = np.zeros(shape =(a,360,720))
    for i in range(360):
        for j in range(720):
            if Region[i,j] > 0:
                mask_boolean[:,i,j] = 0
            else:
                mask_boolean[:,i,j] = 1
    masked_data = ma.masked_array(data, mask = mask_boolean)  
    
    weight = mask_ENDO2D(clat2D, Region)
    a = len(masked_data[:,0,0])
    weighted_data = np.zeros(a)
    for i in range(a):
        weighted_data[i] =(masked_data[i,:,:]*weight).sum()/(weight).sum()
    return weighted_data

def monthly_mean(data):
    '''
    Mean seasonal cycle
    '''
    monthly_mean = np.zeros(12)
    for i in range(12):
        monthly_mean[i] = np.mean(data[i:len(data):12])
    return monthly_mean
