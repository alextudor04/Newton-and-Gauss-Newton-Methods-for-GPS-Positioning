import numpy as np
import matplotlib.pyplot as plt

c = 299792458.0

def solve_newton(satellites, pseudoranges, initial_guess):
    u = initial_guess.copy()
    errors = []
    max_iter = 50
    k = 0

    while k < max_iter:
        F = np.zeros(4)
        J = np.zeros((4,4))

        for i in range(4):
            xi, yi, zi = satellites[i]
            Pi = pseudoranges[i]
            x, y, z, d, = u

            Ri = np.sqrt((x - xi)**2 + (y - yi)**2 + (z - zi)**2)
            F[i] = Ri + c * d - Pi

            J[i, 0] = (x - xi) / Ri
            J[i, 1] = (y - yi) / Ri
            J[i, 2] = (z - zi) / Ri
            J[i, 3] = c

        delta = np.linalg.solve(J, -F)
        u = u + delta
        errors.append(float(np.linalg.norm(delta)))
        k += 1

        if np.linalg.norm(delta) < 1e-8:
            break

    filtered = [(i + 1, e) for i, e in enumerate(errors) if e > 0]
    xs, ys = zip(*filtered)
    plt.figure()
    plt.plot(xs, [np.log10(e) for e in ys], 'bo-', label = 'Newton')
    plt.legend()
    plt.xlabel('Iteration k')
    plt.ylabel('log10(||Δuk||₂)')
    plt.title('Newton Method Convergence')
    plt.grid(True)
    plt.savefig('newton_error_plot.png')
    plt.close()

    return u, len(errors), errors

def solve_gauss_newton(satellites, pseudoranges, initial_guess):
    u = initial_guess.copy()
    errors = []
    N = len(satellites)
    max_iter = 50
    k = 0

    while k < max_iter:
        F = np.zeros(N)
        J = np.zeros((N, 4))

        for i in range(N):
            xi, yi, zi = satellites[i]
            Pi = pseudoranges[i]
            x, y, z, d, = u

            Ri = np.sqrt((x - xi)**2 + (y - yi)**2 + (z - zi)**2)
            F[i] = Ri + c * d - Pi

            J[i, 0] = (x - xi) / Ri
            J[i, 1] = (y - yi) / Ri
            J[i, 2] = (z - zi) / Ri
            J[i, 3] = c

        delta = np.linalg.solve(J.T @ J, -J.T @ F)
        u = u + delta
        errors.append(float(np.linalg.norm(delta)))
        k += 1

        if np.linalg.norm(delta) < 1e-8:
            break

    filtered = [(i + 1, e) for i, e in enumerate(errors) if e > 0]
    xs, ys = zip(*filtered)
    plt.figure()
    plt.plot(xs, [np.log10(e) for e in ys], 'ro-', label = 'Gauss-Newton')
    plt.legend()
    plt.xlabel('Iteration k')
    plt.ylabel('log10(||Δuk||₂)')
    plt.title('Gauss-Newton Method Convergence')
    plt.grid(True)
    plt.savefig('gauss_newton_error_plot.png')
    plt.close()

    return u, len(errors), errors

if __name__ == '__main__':
    import pandas as pd
    df = pd.read_csv('Dataset.csv')
    satellites = df[['SvPositionXEcefMeters',
                     'SvPositionYEcefMeters',
                     'SvPositionZEcefMeters']].to_numpy()

    pseudoranges = df['correctedPrM'].to_numpy()
    satellites_4 = satellites[:4]
    pseudoranges_4 = pseudoranges[:4]

    initial_guess = np.array([0.0, 0.0, 6371000.0, 0.0])

    position_n, iterations_n, errors_n = solve_newton(satellites_4, pseudoranges_4, initial_guess)
    position_gn, iterations_gn, errors_gn = solve_gauss_newton(satellites, pseudoranges, initial_guess)

    # Final Coords Table
    print("=" * 70)
    print(f"{'Method':<20} {'x (m)':>10} {'y (m)':>10} {'z (m)':>10} {'d (s)':>12}")
    print("=" * 70)
    print(f"{'Newton':<20} {position_n[0]:>10.2f} {position_n[1]:>10.2f} {position_n[2]:>10.2f} {position_n[3]:>12.8f}")
    print(f"{'Gauss-Newton':<20} {position_gn[0]:>10.2f} {position_gn[1]:>10.2f} {position_gn[2]:>10.2f} {position_gn[3]:>12.8f}")
    print("=" * 70)

    def calculate_eoc(errors):
        errors = [e for e in errors if e > 0]
        eoc = []
        for k in range(2, len(errors)):
            numerator = np.log(errors[k] / errors[k - 1])
            denominator = np.log(errors[k - 1] / errors[k - 2])
            eoc.append(numerator / denominator)
        return eoc

    # Newton EOC Table
    print("\nNewton Convergence")
    print("=" * 35)
    print(f"{'k':<5} {'||Δuk||':>15} {'EOC':>10}")
    print("=" * 35)
    eoc_n = calculate_eoc(errors_n)
    for k, e in enumerate(errors_n):
        eoc_val = f"{eoc_n[k - 2]:.4f}" if k >= 2 and k - 2 < len(eoc_n) else "N/A"
        print(f"{k + 1:<5} {e:>15.6e} {eoc_val:>10}")
    print("=" * 35)

    # Gauss-Newton EOC Table
    print("\nGauss-Newton Convergence")
    print("=" * 35)
    print(f"{'k':<5} {'||Δuk||':>15} {'EOC':>10}")
    print("=" * 35)
    eoc_gn = calculate_eoc(errors_gn)
    for k, e in enumerate(errors_gn):
        eoc_val = f"{eoc_gn[k - 2]:.4f}" if k >= 2 and k - 2 < len(eoc_gn) else "N/A"
        print(f"{k + 1:<5} {e:>15.6e} {eoc_val:>10}")
    print("=" * 35)

    newton_filtered = [(i + 1, e) for i, e in enumerate(errors_n) if e > 0]
    gauss_filtered = [(i + 1, e) for i, e in enumerate(errors_gn) if e > 0]
    nx, ny = zip(*newton_filtered)
    gx, gy = zip(*gauss_filtered)

    # EOC Graph
    plt.figure()
    plt.plot(nx, [np.log10(e) for e in ny], 'bo-', label='Newton')
    plt.plot(gx, [np.log10(e) for e in gy], 'ro-', label='Gauss-Newton')
    plt.xlabel('Iteration k')
    plt.ylabel('log10(||Δuk||₂)')
    plt.title('Convergence Comparison')
    plt.legend()
    plt.grid(True)
    plt.savefig('combined_plot.png')
    plt.close()
