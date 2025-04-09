from vpython import sphere, vector, rate, label, scene, color
import numpy as np
import glob
import os


def animate_positions(positions_dir="positions_files", fps=100, particle_radius=10.0):
    """
    Animate particle motion from a series of position text files.

    Each file should be formatted with one line per particle: x y z
    Files must be named like 'positions_0000.txt', 'positions_0001.txt', etc.

    Args:
        positions_dir (str): Directory containing position text files
        fps (int): Frames per second for the animation
        particle_radius (float): Radius of the particle spheres
    """

    # --- Configure the scene ---
    scene.title = "Particle Animation"
    scene.width = 1000
    scene.height = 800
    scene.background = color.white
    scene.forward = vector(0, -0.3, -1)
    scene.range = 100  # Adjust if your particles move in a larger domain

    # --- Load and sort all timestep files ---
    files = sorted(glob.glob(os.path.join(positions_dir, "positions_*.txt")))
    if not files:
        raise FileNotFoundError(f"No position files found in '{positions_dir}'")

    # --- Load first file to determine number of particles ---
    positions = np.loadtxt(files[0])
    num_particles = positions.shape[0]

    # --- Create particle spheres ---
    particles = [
        sphere(
            pos=vector(*positions[i]),
            radius=particle_radius,
            color=color.blue,
            make_trail=False,
        )
        for i in range(num_particles)
    ]

    # --- Create label to show current timestep ---
    step_label = label(
        pos=vector(0, 0, 0),
        text='Timestep: 0',
        xoffset=20,
        yoffset=40,
        height=20,
        color=color.black,
        box=False
    )

    # --- Animation loop ---
    for step, filename in enumerate(files):
        data = np.loadtxt(filename)

        if data.shape[0] != num_particles:
            raise ValueError(f"Mismatch in number of particles at timestep {step}.")

        # Update particle positions
        for i in range(num_particles):
            particles[i].pos = vector(*data[i])

        # Update timestep label
        step_label.text = f"Timestep: {step}"

        rate(fps)  # control the animation speed