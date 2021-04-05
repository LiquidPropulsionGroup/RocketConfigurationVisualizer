
from rocketcea.biprop_utils.rho_isp_plot_obj import RhoIspPlot
from rocketcea.biprop_utils.veh_stage_obj import ReferenceStage

rp = RhoIspPlot(bipropL=[('LOX','LH2'), ('LOX','RP1')], Pc=1000., eps=16)


stg_obj = ReferenceStage( WtPayload=1600000.0, volCuInRef=75500.0*1728, WtInertRef=300000.0,
                          Name='Saturn V 1st stg')

rp.add_rho_isp_contours(label_frac_pos=0.2)

rp.add_stage_param_contours( stg_obj, set_param='VolPropellant', param_value=75500.0*1728,
                             label_frac_posD={'GLOW':0.1, 'CubicFt':.4, 'MassFrac':.5},
                             plot_param_valD={'DeltaV':[8000,10000,12000,14000,16000,18000], 
                                              'MassFrac':[0.86,0.88,0.9,0.92,0.94,0.95],
                                              'GLOW':[3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6]    },
                             plot_paramL=['DeltaV','GLOW'], num_ticks=6)

png_name = __file__[:-2] + 'png'
rp.savefig(png_name, dpi=120)
rp.show()

