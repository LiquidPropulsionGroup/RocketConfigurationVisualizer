from rocketprops.rocket_prop import get_prop

# https://rocketprops.readthedocs.io/en/latest/index.html
p = get_prop('rp-1')#c3h8, rp-1
p.plot_sat_props(save_figures=True)