
  #treestring='(((-410.142,1)V,(((T,98.7371)V,(T,p1)V)M,((p0,O)V,(p2,p2)V)M)M)S,(9.23196e-05,0.00913319)V)M'

  #treestring='((-7.55372e-05,0)V,(p2,0)V)M'

  # JUN_SEA_800k
  #treestring= '(((L,Q)V,((T,p1)V,((((((-399.901,p2)V,(-0.000600906,p2)V)M,(p2,0.178955)V)M,(0.514886,0.143141)V)S,(p1,p0)V)S,(p2,p0)V)S)M)A,(-5.86539e-05,4.82462e-05)V)M' # at 80
  #treestring= '(((p0,p0)V,((T,p1)V,((((((0.0121268,p2)V,(15.4748,p2)V)M,(p2,0.163515)V)M,(0.512459,0.0918413)V)S,(p1,p0)V)S,(p2,p0)V)S)M)A,(-5.94995e-05,4.82462e-05)V)M' # at 70

  # simplified
  #treestring= '(((p0,p0)V,((T,p1)V,((((((0.0121268,0.0)V,(15.4748,0.0)V)M,(p2,0.163515)V)M,(0.512459,0.0918413)V)S,(p1,p0)V)S,(p2,p0)V)S)M)A,(-5.94995e-05,4.82462e-05)V)M' # at 70

  # JUN_dust 400

  # SEA LEVEL DATA. SL_t800k_q484_p80_col1_12. At iteration 30
  #treestring = '((0.00686355,0.000170046)V,((0.00785789,p0)V,(((p1,-0.0745519)V,(0.910784,-0.766288)V)A,((p2,T)V,(-0.000190906,p1)V)S)A)M)M'

  # SEA LEVEL FINAL ITERATION
  #treestring = '((-9.04601e-05,-1.51439e-05)V,(((T,p0)V,((p1,-0.876754)V,(p3,4.97253)V)S)M,((-0.00729842,p0)V,((58.2551,p1)V,((T,T)V,((p1,T)V,((p3,T)V,((-27.4471,p0)V,((p2,T)V,(T,T)V)S)A)M)A)A)S)M)A)M'

  treestring = '((7.20581e-05,7.08437e-05)V,((((p2,p0)V,((0.00230651,-0.00859469)V,((p1,-0.0916722)V,((p2,T)V,(T,T)V)A)A)M)S,(((p1,T)V,((p0,p0)V,(-0.0432388,T)V)M)A,(T,p1)V)S)S,((1,p0)V,(p0,T)V)M)S)M'

  # EPICA. EP_t800k_q484_p80_col1_12. COL 1 At iteration 30
  #treestring = '(((p1,0.00450953)V,((T,-0.0176874)V,(((-0.102754,0.431628)V,(0.00521597,T)V)S,((-0.00489772,p0)V,((p0,p1)V,(p3,T)V)S)A)M)M)M,(p1,-0.149613)V)M'

  # EPICA pow4 cost function. Latest.
  #treestring = '((0.000115338,-9.78074e-06)V,((((T,T)V,((p3,p0)V,(T,p0)V)S)M,(1,p0)V)A,(((p3,T)V,((p1,p1)V,((2.86446,T)V,((T,T)V,((207.884,p1)V,(T,p0)V)M)M)A)M)S,((p0,1.23357)V,((T,T)V,((((((p2,T)V,(T,p1)V)M,((p1,T)V,(((0,p0)V,(T,p1)V)S,(p3,p1)V)M)A)S,(p3,p1)V)A,((T,p1)V,(T,T)V)A)M,(((p1,T)V,(T,T)V)S,(T,p1)V)M)S)M)S)S)M)M'

  # TEST TREES
  #treestring= '(((0,p1)V,(0,0)V)A,((0,p0)V,(0,1)V)A)M'
  #treestring= '((0,p0)V,((0,p0)V,(0,1)V)A)M'

  # this becomes p1**2. it seems that the right multiplicant is copied into the left
  #treestring= '(((0,p0)V,(0,1)V)A,((0,p1)V,(0,0.5)V)A)M'

  # this works and gives 0:
  #treestring= '(((0,0)V,(0,0)V)A,(0,p1)V)M'

  # when both top multiplicants contain 2-arg children, right top arg is copied into left

  #treestring= '((0,p1)V,(0,p1)V)M'

