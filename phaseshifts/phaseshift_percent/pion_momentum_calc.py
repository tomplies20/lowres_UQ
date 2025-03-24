import numpy as np

hbarc = 197.326
M = (938.272 + 939.565) / 2.0  # averaged neutron/proton mass in MeV
units_factor = hbarc * hbarc / M
def Elab(p):
    return 2 * p ** 2 * hbarc ** 2 / M

print(Elab(139/hbarc))
print(139/hbarc/(600/hbarc))
print(Elab(1/hbarc))
print(Elab(2.8))

def mom(E):
    return np.sqrt(M * E / 2 / hbarc ** 2)

print(mom(394.6539576))

print(2*hbarc)