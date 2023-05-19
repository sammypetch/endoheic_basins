
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

