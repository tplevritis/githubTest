import numpy as np
import matplotlib.pyplot as plt

def panel_method(x_coords, y_coords, alpha, V):
    num_panels = len(x_coords) - 1

    panel_length = np.sqrt(np.diff(x_coords)**2 + np.diff(y_coords)**2)
    panel_angle = np.arctan2(np.diff(y_coords), np.diff(x_coords))
    panel_center_x = 0.5 * (x_coords[:-1] + x_coords[1:])
    panel_center_y = 0.5 * (y_coords[:-1] + y_coords[1:])

    A = np.zeros((num_panels, num_panels))
    B = np.zeros(num_panels)

    for i in range(num_panels):
        for j in range(num_panels):
            if i == j:
                A[i, j] = 0.5 * np.pi
            else:
                delta_x = panel_center_x[i] - x_coords[j + 1]
                delta_y = panel_center_y[i] - y_coords[j + 1]
                r_squared = delta_x**2 + delta_y**2
                r = np.sqrt(r_squared)
                cos_theta = delta_x / r
                sin_theta = delta_y / r
                cos_phi = np.cos(panel_angle[i] - alpha)
                sin_phi = np.sin(panel_angle[i] - alpha)

                A[i, j] = cos_theta * cos_phi + sin_theta * sin_phi
                B[i] -= -2 * np.pi * V * sin_phi * panel_length[j]

    gamma = np.linalg.solve(A, B)
    Cl = 2 * np.sum(gamma * panel_length) / V
    Cd = 2 * np.sum(gamma * panel_length * np.sin(panel_angle - alpha)) / V

    return Cl, Cd

# Airfoil coordinates (Example: NACA 2412 airfoil)
x_coords = np.linspace(0, 1, 101)
y_coords = 0.12 * (0.2969 * np.sqrt(x_coords) - 0.1260 * x_coords - 0.3516 * x_coords**2 + 0.2843 * x_coords**3 - 0.1015 * x_coords**4)

# Air properties
alpha = np.radians(5)  # Angle of attack in radians
V = 50.0  # Free-stream velocity (m/s)

# Calculate coefficients using panel method
Cl, Cd = panel_method(x_coords, y_coords, alpha, V)

# Output results
#print(f"Angle of Attack (α): {np.degrees(alpha)} degrees")
print(f"Lift Coefficient (Cl): {Cl:.4f}")
print(f"Drag Coefficient (Cd): {Cd:.4f}")

# Plot the airfoil shape
plt.figure(figsize=(8, 4))
plt.plot(x_coords, y_coords, 'b-')
plt.gca().invert_yaxis()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Airfoil Shape')
plt.grid()
plt.show()
