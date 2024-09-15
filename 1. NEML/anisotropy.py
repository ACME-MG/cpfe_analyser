
# YOUNGS   = 151875
# POISSONS = 0.365
# SHEAR    = 61500

# C_1 = YOUNGS*(1-POISSONS)/(1+POISSONS)/(1-2*POISSONS)
# C_2 = YOUNGS*POISSONS/(1+POISSONS)/(1-2*POISSONS)
# C_3 = 2*SHEAR

# print(f"C_1 = {C_1}")
# print(f"C_2 = {C_2}")
# print(f"C_3 = {C_3}")

C_1 = 250000
C_2 = 151000
C_3 = 123000

print(f"A = {2*C_3/(C_1-C_2)}")
