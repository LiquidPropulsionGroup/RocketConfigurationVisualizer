from rocketcea.cea_obj import CEA_Obj

C = CEA_Obj(oxName='LOX', fuelName='RP1')

full_cea = C.get_full_cea_output(Pc=15.2, MR=2.3, eps=None, PcOvPe=15.2/0.9, pc_units='bar', output='KJ', short_output=1 )
print('this is the cea string:------------------------------------------------------------------------------------------------------------------')
print(full_cea)
