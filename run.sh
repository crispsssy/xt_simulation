#gasFile=gasfile/He_90.0000_iC4H10_10.0000_H2O_0.0000_1.000000atmB0.gas
gasFile=gasfile/He_90_iC4H10_10_760Torr_293.15K_10emin_to_500000emax.gas
ionMob=ionMobility/IonMobility_He+_He.txt
geoFile=geofile/prototypeGeo_9anode.root
HV=1800
evenum=100
fileId=002
prim_particle=muon
prim_momentum=1e9


build/xt_sim ${gasFile} ${ionMob} ${geoFile} ${HV} ${evenum} ${fileId} ${prim_particle} ${prim_momentum}
