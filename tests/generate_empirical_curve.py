import numpy as np

bp_points = [
    10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000,
    900, 800, 700, 600, 500, 400, 300, 200, 100, 50, 25
]
cm_points = [
    1.20, 1.45, 1.75, 1.95, 2.25, 2.65, 3.25, 3.75, 4.15, 4.50,
    4.65, 4.85, 5.05, 5.30, 5.60, 6.00, 6.50, 7.30, 8.50, 9.30, 9.60
]

log_bp = [np.log10(bp) for bp in bp_points]
# Reversing so x is strictly increasing
x_data = log_bp[::-1]
y_data = cm_points[::-1]

test_bps = [10000, 3000, 2000, 1000, 500, 200, 100, 50]
for bp in test_bps:
    calc_cm = np.interp(np.log10(bp), x_data, y_data)
    print(f"{bp:5d} bp: {calc_cm:.2f} cm")
