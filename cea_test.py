from rocketcea.cea_obj import CEA_Obj

C = CEA_Obj(oxName='LOX', fuelName='RP1')

full_cea = C.get_full_cea_output(Pc=100, MR=1.0, eps=40, subar=None, PcOvPe=None, show_transport=1, pc_units='bar', output='KJ', short_output=1 )
print('this is the cea string:------------------------------------------------------------------------------------------------------------------')
print(full_cea)
