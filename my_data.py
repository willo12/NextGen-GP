import os
import numpy as np
from scipy.interpolate import interp1d
import warnings

#import random

def norm(series):
#  return series

  if (len(series) == 0):
    raise Exception("Error: time series must be non-empty")

  new_series = series - np.mean(series)   

  return 2*(new_series)/max(abs(new_series))

def norm_lambda(series):
#  return series

  return lambda S: 2*(S - np.mean(series))/max(abs(series - np.mean(series)))



def norm_error(series):
    """Normalize time series.
    """
#  return series

    new_series = deepcopy(series)
    new_series[:,0] = series[:,0] - np.mean(series[:,0])   

    return 2*(new_series)/max(abs(new_series[:,0]))

def get_grid(obs_times,dt = 50,t0 = -int(2.5e6),t_end=0):

  obs_times = np.round(obs_times).astype(int)

  runsteps = int((t_end - t0)/dt+1)

  indices=[]
  for step in xrange(runsteps):

    t_synth = int(round(t0 + step*dt))
    
    I=np.where(obs_times==t_synth)[0]  
    if len(I)>0:     # detect a common time point
     
      indices.append(step)
    else:
#      print "missing: step=%d for t_synth=%d \n"%(step,t_synth)
      pass

#  print 'Prepared time indices of length %d vs obs lengths %d'%(len(indices),len(obs_times))

  return np.array(indices)

def dlmread(filepath, datatype='float'):
  """
  Read space delimited data. Select 'float' or 'int' for datatype.
  """
  str_data = []
  data = []

#  print filepath

  fobj = open(filepath,'r')
  str_data = fobj.readlines()
  fobj.close()

  if datatype == 'float':
    length = 0
    for eachline in str_data:
      for elem in eachline.split():
        data.append(float(elem))  
        length=length+1
  elif datatype == 'int':
    length = 0
    for eachline in str_data:
      for elem in eachline.split():
        data.append(int(elem))  
        length=length+1
  else:
    length = 0
    for eachline in str_data:
      for elem in eachline.split():
        data.append(elem)  
        length=length+1

  width = len(eachline.split())
  data = np.array(data).reshape(np.floor(length/ width),width)
  
  return data



def get_forcing(forcing_dir=os.environ['HOME'] + '/PROJECTS/paper_glac/',forcing_files = 'j_65north_trunc.txt^2', random_comp = 0.0 ):
  """
  Reads forcing file, normalizes forcing data and interprets positive or negative time axis.

  Time axis assumed to be column 0 of first forcing file.

  example: forcing_files 'j_65north_trunc.txt*1e3^2-3_and_milankovitch.data.txt^1-2' picks cols 2 and 3 from first file, and 1,2 from second. It also interprets time in kyr

  """

  if random_comp>0:
    import random

 # forcing files to format: [('milankovitch.data.txt',(1,2)),]  

  filepluscols = forcing_files.split('_and_')
  
  forcing_files = []
  for fpc in filepluscols:
    fpc_split = fpc.split('^')

    fname_part = fpc_split[0]
    L = fname_part.split('*')
    fname_part = L[0]
    if len(L)>1:
      time_fact = int(round(float(L[1])))
    else:
      time_fact = 1

    cols_str = fpc_split[1]

    cols = tuple([int (e) for e in cols_str.split('-')])
    forcing_files.append((fname_part, cols))


#  raise Exception('ok, i received forcing files: %s'%forcing_files)

  ffs = [] # forcings

  first_flag=True
  for fpair in forcing_files:
    cols=fpair[1] # columns
    ff=fpair[0] # file

    tmp = np.loadtxt(forcing_dir + ff )
    if first_flag: # set time axis once, for the first forcing file only: so first forcing file time axis is used
      if tmp[1,0] >0: # time axis positive in past direction
        T = -tmp[:,0]*time_fact
      else:
        T = tmp[:,0]*time_fact

      ffs.append(T)
      first_flag=False
     
    if random_comp>0:

      for i in cols:
        ffs.append( tmp[:,i] + random.uniform(-random_comp, random_comp) )

    else:

      for i in cols:
        ffs.append( tmp[:,i]  )


  del tmp 
  ffs = zip(*ffs)
  ffs.sort()

 # print("max yo %g"%np.max(np.array(ffs)[1,:]))

  return np.array(ffs)



def get_files(forcing_files = 'j_65north_trunc.txt^2' , obs_file ='PROJECTS/paper_glac/raymo_d18_trunc.txt',obs_error=-1,detrend=False,cols=[5,], home = os.environ['HOME'], t_start = -int(2.5e6),forcing_substeps=2, random_comp = 0.0):
  """
  Example: obs_file = 'PROJECTS/paper_glac/epica_T.txt__mul__-1' multiplies series by -1
  """

  L = obs_file.split("__mul__")
  obs_file = L[0]
  if len(L)>1:
    mul_fac = float(L[1])
  else:
    mul_fac = 1.0
	
  raw_obs=np.loadtxt(os.path.join(home,obs_file))

  raw_obs[:,1] = mul_fac*raw_obs[:,1]

  if (abs(raw_obs[2,0] - raw_obs[1,0]) <= 10):   # REMOVE THIS LATER, AND REQUIRE FILE STANDARD IN YRS 
    raw_obs[:,0] = raw_obs[:,0]*1000 # from kyr to yr
    print("Data not in yrs, ts=%.2f, converting from assumed kyr."%abs(raw_obs[2,0] - raw_obs[1,0]))
#  else:
#    otime = raw_obs[:,0] # it's already in yr

  if (raw_obs[2,0] > 0):
    raw_obs[:,0] = -raw_obs[:,0] # ensure negative time axis values, will use sort later
    raw_obs = np.flipud(raw_obs)

  obs = raw_obs

#  if detrend is True:
#    model = np.polyfit(otime,d18o, 2)
#    fit = np.polyval(model, otime)

#    d18o = d18o - fit
#    print 'Detrended obs.'

#  obs = zip(otime,d18o,dust)
 
#  obs.sort() 
 
#  obs = np.array(obs)

  I_obs_start = np.where(obs[:,0]==t_start)  

  if (len(I_obs_start)>0):
    i_obs_start = I_obs_start[0]
    if not isinstance(i_obs_start,int):
        if len(i_obs_start) > 0:
          i_obs_start = i_obs_start[0]
          print("t_start=%d found in obs at index=%d"%(t_start, i_obs_start))
        else:
          i_obs_start = 0
          warnings.warn("Warning! t_start=%d not found in obs, setting i_obs_start=%d"%(t_start,i_obs_start))
    else:
      print("t_start=%d found in obs at index=%d"%(t_start, i_obs_start))

  else:
    i_obs_start = 0
    warnings.warn("Warning! t_start=%d not found in obs, setting i_obs_start=%d"%(t_start, i_obs_start))

  obs = obs[i_obs_start:,:]

  if obs_error > 0:
    obs[:,1] = norm(obs[:,1]-obs_error*obs[:,2])  
    obs[:,2] = norm(obs[:,1]+obs_error*obs[:,2])  

  else:
    for i in range(1,obs.shape[1]):
      obs[:,i] = norm(obs[:,i])

    if obs_error == -1:
      obs = obs[:,:-1]

  ffs = get_forcing(forcing_files = forcing_files, random_comp = random_comp)  # expecting data in yr

#  raise Exception("yo %d"%ffs.shape[1])

  # interpolate forcing

#  print ffs[0]

  I_forc_start = np.where(ffs[:,0]==t_start)  

  if len(I_forc_start)>0:

    try:
        i_forc_start = I_forc_start[0]
        if not isinstance(i_forc_start,int):
            if len(i_forc_start)>0:
              i_forc_start = i_forc_start[0]        
              print("t_start=%d found in forcing at index=%d"%(t_start, i_forc_start))
            else:
              i_forc_start = 0
              # should we raise an exception here?
              warnings.warn("Warning! t_start=%d not found in forcing, setting index=%d. Check length forcing file. WRONG LENGTH COULD YIELD A SEG FAULT!"%(t_start, i_forc_start) )


    except:
        raise Exception("Problem at ffs: %f"%ffs[0,0])

  else:
    i_forc_start = 0
    warnings.warn("Warning! t_start=%d not found in forcing, setting index=%d. Check length forcing file."%(t_start, i_forc_start))

  slice_ffs_time = slice(i_forc_start,None)  # single slice for time dimension

  slices_ffs =  ( slice_ffs_time ,  slice(None)   )  # tuple containing all slices
#  slices_ffs[0] = slice_ffs_time

#  raise Exception(slices_ffs)


  ffs = ffs[ slices_ffs ]

  dt_forc_raw = abs(ffs[1,0] - ffs[0,0])
  int_step = int(round(dt_forc_raw/forcing_substeps))

  t_end = int(round(ffs[-1,0]))

  # finer times to interpolate to
  itime = np.arange(int(round(ffs[0,0])),t_end+int_step ,int_step )

#  print itime

  sh = ffs.shape

  Iffs = np.zeros((len(itime),sh[1]))

  Iffs[:,0] = itime

  ftime = np.squeeze(ffs[:,0]) 
  for iforc in xrange(1,sh[1]):
     
    Iffs[:,iforc] = interp1d(ftime , np.squeeze(ffs[:,iforc]) )(itime)

  for i in range(1,Iffs.shape[1]):  # rescale forcing data to between -2 and 2
    Iffs[:,i] = norm(Iffs[:,i])

  dt_forc = abs(Iffs[1,0] - Iffs[0,0])

  I = get_grid(np.squeeze(obs[:,0]), dt=dt_forc, t0 = t_start, t_end = t_end) # I subset of Iffs times, must have same length as obs[:,0]

  if I.shape[0] > obs.shape[0]: # c code will compare result[I[i]] values to obs[i], so must match
    raise Exception("Data Setting Error: I must be shorter than obs. All obs time steps must be divisible by dt_forc. Also, t_start (%d vs obs %d) must be within data file range. Shape I: %s, shape obs: %s."%(t_start, obs[0,0], I.shape[0] , obs.shape[0]))

  obs = obs[:I.shape[0],:]

  return (Iffs, obs,I)

def glacial(smooth=0,t_obs=0,t_forc=0, t_start = -int(2.5e6), dt=50, startscore_t=None, ts_factor=2, interp=False, obs_file ='PROJECTS/paper_glac/raymo_d18_trunc.txt', forcing_files = 'j_65north_trunc.txt*1e3^2',forcing_substeps=2, random_comp=0.0):


  Iffs, obs_tmp,I = get_files(forcing_files = forcing_files, obs_file = obs_file, t_start = t_start,forcing_substeps=forcing_substeps, random_comp = random_comp)

  Iffs = Iffs[t_forc:,:]
#  I = get_grid(obs_tmp[:,0],dt,t0 = t_start)


  obs_tmp = obs_tmp[t_obs:,:] # truncate time

  obs0 = obs_tmp.shape[0]
  obs1 = obs_tmp.shape[1]

  time_obs = obs_tmp[:,0] # time of obs
  ts = obs_tmp[:,1]  # the actual time series

  if smooth>0:
    ts_tmp = smoothGaussian(obs_tmp[:,1],degree=smooth)
    offset = (len(ts) - len(ts_tmp))/2
    ts[offset:len(ts)-offset-1] = ts_tmp  

  dts = np.zeros(len(ts));dts[:-1]=ts[1:]-ts[:-1]

  if interp is True:
    f = interp1d(time_obs, ts)
    time_fine = np.arange(t_start,dt,dt)
    ts_fine = f(time_fine)

    f = interp1d(time_obs, dts)
    dts_fine = f(time_fine)

    obs0 = len(ts_fine)
    obs1 = 3

    obs = np.zeros((obs0,obs1))
#  obs[:,0] = obs_tmp[:,0]
#  obs[:,1] = ts
#  obs[:,2] = dts
##  obs[:,3] = obs_tmp[:,2]

    obs[:,0] = time_fine
    obs[:,1] = ts_fine
    obs[:,2] = dts_fine

  else:
#    obs0 = obs_tmp.shape[0]
#    obs1 = obs_tmp.shape[1]

#    obs = np.zeros((obs0,obs1))

#    obs[:,0] = obs_tmp[:,0]
#    obs[:,1] = ts
#    obs[:,2] = dts

#  I = get_grid(obs[:,0],dt,t0 = t_start) # overlap with obs, should be all indices now

    obs = obs_tmp

  if startscore_t is not None:

    startscore_i=np.where(obs[:,0]==startscore_t)[0]  
    if len(startscore_i)>0:     # detect a common time point
      startscore_i = startscore_i[0]
    else:

      if (startscore_t < obs[0,0]):
        warnings.warn("Warning: startscore_t (value %d) < earliest obs obs[0,0] (value %d). Setting startscore_i=0"%(startscore_t, obs[0,0]))
      else:
        # should we raise an exception here?
        warnings.warn("Not found: startscore_i, even though startscore_t (value %d) > earliest obs obs[0,0] (value %d). Check obs file ok. Setting startscore_i=0. Answer will likely be wrong, and using wrong length files can cause seg faults.")

      startscore_i = 0

  else:
    startscore_i=0

  return Iffs, obs,I, startscore_i

def epica(forcing_files = [('j_65north_trunc.txt', (1,))  ], obs_file ='PROJECTS/paper_glac/epica_T.txt',detrend=False,cols=[5,], dt_forc=50, home = os.environ['HOME'], ts_factor=1, startscore_t=None,smooth=0):
	
  obs=np.loadtxt(os.path.join(home,obs_file))

  otime = obs[:,0]
  d18o = obs[:,1]

  if detrend is True:
    model = np.polyfit(otime,d18o, 2)
    fit = np.polyval(model, otime)

    d18o = d18o - fit
    print 'Detrended obs.'

  obs[:,1] = norm(obs[:,1])

  if smooth>0:
 
    ts_tmp = smoothGaussian(obs[:,1],degree=smooth)

    offset = (len(obs[:,1]) - len(ts_tmp))/2

    obs[offset:len(obs[:,1])-offset-1,1] = ts_tmp  


  ffs = get_forcing(forcing_files = forcing_files)

  # interpolate forcing

#  print ffs[0]

  itime = np.arange(-2.5e3,0.1,0.1)

#  itime = np.arange(obs[0,0],dt_forc,dt_forc)*1e-3  # time values to interpolate to. Forcing is in 1000 yrs

  sh = ffs.shape

  Iffs = np.zeros((len(itime),sh[1]))

  Iffs[:,0] = itime*1e3

  ftime = np.squeeze(ffs[:,0]) 
  for iforc in xrange(1,sh[1]):
     
    Iffs[:,iforc] = interp1d(ftime , np.squeeze(ffs[:,iforc]) )(itime)

  I = get_grid(np.squeeze(obs[:,0]) ,dt = 50,t0 = obs[0,0],t_end=0 )

#  return (Iffs, obs,I, ts_factor)

  if startscore_t is not None:
    startscore_i=np.where(obs[:,0]==startscore_t)[0]  
    if len(startscore_i)>0:     # detect a common time point
      startscore_i = startscore_i[0]
    else:
#      print "Not found"
      pass

  else:
    startscore_i=0

  return Iffs, obs,I, ts_factor, startscore_i


def onseries(factor=4,smooth=0):

  T = np.loadtxt('T')+1;  
  n = T.shape[0]

  tim = np.zeros((n*factor,1))
  tim[:,0] = np.arange(0,n*factor)
  t_obs = np.arange(0,n*factor,factor)

  obs = np.zeros((n,3))
  obs[:,0] = t_obs
  obs[:,1] = T[:,0]

  ts = obs[:,1]
  dts = np.zeros(len(ts))
  dts[:-1] = ts[1:] - ts[:-1]
  dts_smt_trunc=smoothGaussian(dts,degree=smooth)
  dts_smt=dts;

  offset = (len(dts_smt) - len(dts_smt_trunc))/2

  dts_smt[offset:len(dts_smt)-offset-1] = dts_smt_trunc;

  obs[:,2] = dts_smt
  ts_factor = 1

  return tim, obs, t_obs, ts_factor

def synth(forcing_files = [('j_65north_trunc.txt', (1,))  ], obs_file ='PROJECTS/paper_glac/synth.txt',cols=[5,], dt_forc=50, home = os.environ['HOME'], ts_factor=1):
	
  obs=np.loadtxt(os.path.join(home,obs_file))

  otime = obs[:,0]
  d18o = obs[:,1]

#  obs[:,1] = norm(obs[:,1])

  ffs = get_forcing(forcing_files = forcing_files)

  # interpolate forcing

#  print ffs[0]

  itime = np.arange(obs[0,0],dt_forc,dt_forc)*1e-3  # time values to interpolate to. Forcing is in 1000 yrs

  sh = ffs.shape

  Iffs = np.zeros((len(itime),sh[1]))

  Iffs[:,0] = itime*1e3

  ftime = np.squeeze(ffs[:,0]) 
  for iforc in xrange(1,sh[1]):
     
    Iffs[:,iforc] = interp1d(ftime , np.squeeze(ffs[:,iforc]) )(itime)

  I = get_grid(np.squeeze(obs[:,0]) ,dt = 50,t0 = obs[0,0],t_end=0 )

  return (Iffs, obs,I, ts_factor)



