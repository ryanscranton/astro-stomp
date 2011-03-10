n_z = 49


;load_dp,/x
file_output = 'STOMP_OUTPUT/output_photoqso_low_0_1Mpc.fit'
;file_output = 'STOMP_OUTPUT/output_LRG_MgII.fit'

;seed_path = '/mnt/'
;X_file = seed_path + 'raid-cita/menard/DATA/SDSS/LRG/megazlrg.fit'
;target_file = seed_path+'raid-cita/menard/DATA/SDSS/Quasar/Spectro/dr7qso.fit'
;target_file = 'QSO_MgII_Nestor.fit'


;X = MRDFITS(X_file,1)
;target = MRDFITS(target_file,1)
stomp = MRDFITS(file_output,1)


z_min = 0. ;min(stomp.z_target)
z_max = 3.0 ;max(stomp.z_Target)
z_vector = make_vector_structure(z_min,z_max,n_z)
neighbor = REPLICATE( {mean:0.d0, mean_err:0.d0, median:0.d0}, n_z)

!p.multi=[0,3,3]

for i=0,n_z-1 do begin
   index_in_bin = where( (stomp.z_target gt z_vector.bound(i)) AND $
                         (stomp.z_target lt z_vector.bound(i+1)) )
   if (index_in_bin)(0) gt 0 then begin
      neighbor(i).mean = AVG(stomp(index_in_bin).n_neighbor)
      neighbor(i).mean_err = STDDEV((stomp(index_in_bin).n_neighbor))/sqrt(1.*n_elements(index_in_bin))
      neighbor(i).median = MEDIAN(stomp(index_in_bin).n_neighbor)
;   plothist,stomp(index_in_bin).n_neighbor
   endif

endfor


normalization = ANGDIST_LAMBDA(z_vector.mean,0.,h=0.7)

erase
multiplot,/default
multiplot,[1,2]

y = neighbor.mean * normalization^2
y_err = neighbor.mean_err * normalization^2
y = y - MEDIAN(y)
my_xr = minmax(z_vector.bound)
my_yr = 1.2*max(y)*[-1,1]

plot,z_vector.mean,y,xr=my_xr,ystyle=1,yr=my_yr,$
   ytit='recovered N(z)'
oploterr,z_vector.mean,y,y_err
oplot,my_xr,[0,0],linestyle=2;,color=getcolor('gray',2)
multiplot

;my_bin = 0.01
;my_color=getcolor('black',1)
;PLOTHIST,QSO_phot.redshift,bin=my_bin,peak=1,xr=my_xr,color=my_color,/fill,fcolor=my_color
;PLOTHIST,QSO_phot.z,bin=my_bin,peak=1,xr=my_xr,color=my_color,/fill,fcolor=my_color

;PLOTHIST,target.redshift_abs,bin=my_bin,/fill,/fline,peak=1,xr=my_xr;,/overplot

;PLOTHIST,X.z_phot,bin=my_bin,peak=1,xr=my_xr,color=my_color,/fill,fcolor=my_color,ytit='photo z',xtit='z'

multiplot,/default
forprint, z_vector.mean, y, y_err,textout='output_photoqso_low.out'

end
