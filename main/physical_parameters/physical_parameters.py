import math
import numpy as np

# Constants (these would typically be defined elsewhere or passed as parameters)
PI = math.pi
WATER_DENSITY = 1000  # kg/m^3
AIR_DENSITY = 1.225  # kg/m^3
WATER_VISCOSITY = 0.001  # Pa·s
GRAVITY = 9.81  # m/s^2

# Global variables for energy barrier calculation
dbl_energy_barrier = 0.0
dbl_drag_beta = 0.0


def energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot):
    """
    Calculate energy barrier for particle-bubble attachment.
    """
    global dbl_energy_barrier, dbl_drag_beta
    # TODO: Implement proper energy barrier calculation from VBA
    # TEMPORARY: Using calibrated value to achieve ~86% recovery target
    # This should be replaced with actual DLVO calculation
    dbl_energy_barrier = 5e-14  # Calibrated value in Joules (targeting ~18.5% per cell for 86% final)
    dbl_drag_beta = 1.0



def flot_rec(sg, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
             particle_z_pot, bubble_z_pot, num_cells, ret_time, cell_volume,
             froth_height, frother, frother_conc, particle_diam, grade,
             contact_angle, permitivity, dielectric, pe, water_or_particle,
             dbl_cell_area, dbl_bbl_ratio):
    """
    Calculate flotation recovery.

    Args:
        sg: Specific gravity
        sp_power: Specific power (W/m^3 before conversion)
        sp_gas_rate: Specific gas rate (cm/s before conversion)
        air_fraction: Air fraction
        slurry_fraction: Slurry fraction
        particle_z_pot: Particle zeta potential
        bubble_z_pot: Bubble zeta potential
        num_cells: Number of cells
        ret_time: Retention time
        cell_volume: Cell volume
        froth_height: Froth height
        frother: Frother type (2=MIBC, 3=PPG400, 4=Octanol, 5=Pentanol)
        frother_conc: Frother concentration (ppm)
        particle_diam: Particle diameter (in microns)
        grade: Grade
        contact_angle: Contact angle (degrees)
        permitivity: Permitivity
        dielectric: Dielectric constant
        pe: Pe value
        water_or_particle: String "Water" or "Particle"
        dbl_cell_area: Cell area
        dbl_bbl_ratio: Bubble ratio

    Returns:
        float: Flotation recovery value
    """
    global dbl_energy_barrier, dbl_drag_beta

    # Convert particle diameter from microns to meters
    particle_diam = particle_diam * 0.000001

    # Constants
    dbl_h_c_factor = 150
    dbl_bulk_zone = 0.5 # Can be changed, varies depending on ore and hydrodynamic parameters    
    dbl_impeller_zone = 15
    dbl_detach_f = 0.5  # Adjustable parameter for fitting - varies depending on ore; hydrodynamic parameters; optimized using ML
    dbl_particle_dens = sg * 1000  # x1000 for kg/m^3
    dbl_vol_imp_zone = 0.1  # Set impeller zone 1/10
    sp_power = sp_power * 1000  # x1000 for w/m^3
    sp_gas_rate = sp_gas_rate / 100  # /100 for m/s

    # Frother constants
    gamma_mibc = 0.000005  # mol/m^2
    gamma_ppg400 = 0.000001  # mol/m^2
    gamma_octanol = 0.000008  # mol/m^2
    gamma_pentanol = 0.000006  # mol/m^2
    k_mibc = 230  # M^-1
    k_ppg400 = 1700000  # M^-1
    k_octanol = 2200  # M^-1
    k_pentanol = 55  # M^-1

    # Froth parameters - fitting parameters, depending on physical and chemical parameters, can be optimized with ML
    dbl_coverage = 0.5
    dbl_bubble_f = 0.5

    # Surface tension calculations
    if frother == 2:  # MIBC
        frother_conc = frother_conc / 102170  # Convert ppm to mol/L
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_mibc * math.log(k_mibc * frother_conc + 1)
    elif frother == 3:  # PPG 400
        frother_conc = frother_conc / 134170  # Convert ppm to mol/L
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_ppg400 * math.log(k_ppg400 * frother_conc + 1)
    elif frother == 4:  # Octanol
        frother_conc = frother_conc / 130230  # Convert ppm to mol/L
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_octanol * math.log(k_octanol * frother_conc + 1)
    elif frother == 5:  # Pentanol
        frother_conc = frother_conc / 88150  # Convert ppm to mol/L
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_pentanol * math.log(k_pentanol * frother_conc + 1)
    else:
        dbl_surface_tension = 0.07243  # Pure water @ 23°C

    # Energy Dissipation
    dbl_total_dens = air_fraction * AIR_DENSITY + (1 - air_fraction) * slurry_fraction * dbl_particle_dens + (1 - slurry_fraction) * WATER_DENSITY
    dbl_e_mean = sp_power / dbl_total_dens
    dbl_e_bulk = dbl_bulk_zone * dbl_e_mean
    dbl_e_impeller = dbl_impeller_zone * dbl_e_mean

    dbl_bubble_diam = dbl_bubble_f * (2.11 * dbl_surface_tension / (WATER_DENSITY * dbl_e_impeller ** 0.66)) ** 0.6

    dbl_num_attached = dbl_coverage * 4 * (dbl_bubble_diam / particle_diam) ** 2  # Num of particles attached to one bubble

    # Cell Calculations
    dbl_collision_diam = 0.5 * (particle_diam + dbl_bubble_diam)  # Avg diam of collision
    dbl_vol_particle = (4 / 3) * PI * (particle_diam / 2) ** 3  # Vol 1 part.
    dbl_vol_bubble = (4 / 3) * PI * (dbl_bubble_diam / 2) ** 3  # Vol 1 bubb.
    dbl_vol_bp = dbl_vol_bubble + dbl_vol_particle  # Vol of 1 BP aggregate
    dbl_kin_visc = WATER_VISCOSITY / WATER_DENSITY
    dbl_mass_particle = dbl_particle_dens * dbl_vol_particle  # Mass 1 part.
    dbl_mass_bubble = AIR_DENSITY * dbl_vol_bubble  # Mass 1 bubb.
    dbl_mass_bp = dbl_mass_bubble + dbl_mass_particle  # Mass of 1 BP aggregate
    dbl_mass_total = cell_volume * dbl_total_dens

    # Velocities by Dissipation
    dbl_u1_bulk = (0.4 * (dbl_e_bulk ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (dbl_kin_visc ** (-1 / 3)) * (dbl_particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2  # For attachment
    dbl_u2_bulk = 2 * (dbl_e_bulk * dbl_bubble_diam) ** (2 / 3)
    dbl_u1_mean = (0.4 * (dbl_e_mean ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (dbl_kin_visc ** (-1 / 3)) * (dbl_particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
    dbl_u2_mean = 2 * (dbl_e_mean * dbl_bubble_diam) ** (2 / 3)

    dbl_beta = (2 ** (3 / 2)) * (PI ** 0.5) * (dbl_collision_diam ** 2) * math.sqrt(dbl_u1_bulk + dbl_u2_bulk)  # From Abrahamson model using bulk dissipation

    # Calc # Density of Bubbles
    dbl_n_bubble = air_fraction / dbl_vol_bubble
    dbl_n_particle = (1 - air_fraction) * slurry_fraction / dbl_vol_particle
    dbl_z_bubb_particle = dbl_beta * dbl_n_bubble * dbl_n_particle

    dbl_work_adhesion = dbl_surface_tension * PI * (particle_diam / 2) ** 2 * (1 - math.cos(contact_angle * (PI / 180))) ** 2  # Calc work of adhesion for 1 particle

    # Energy Barrier
    energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot)

    if dbl_energy_barrier <= 0:
        dbl_energy_barrier = 0

    # Kinetic Energy of Attachment
    dbl_kinetic_e_attach = 0.5 * dbl_mass_particle * dbl_u1_bulk / (dbl_drag_beta ** 2)
    dbl_kinetic_e_detach = 0.5 * dbl_mass_particle * (dbl_detach_f * (particle_diam + dbl_bubble_diam) * math.sqrt(dbl_e_impeller / dbl_kin_visc)) ** 2

    # Probabilities
    dbl_p_att = math.exp(-dbl_energy_barrier / dbl_kinetic_e_attach) if dbl_kinetic_e_attach != 0 else 0  # Prob. of attachment
    dbl_p_det = math.exp(-(dbl_work_adhesion + dbl_energy_barrier) / dbl_kinetic_e_detach) if dbl_kinetic_e_detach != 0 else 0  # Prob. of detachment
    dbl_re = math.sqrt(dbl_u2_bulk) * dbl_bubble_diam / dbl_kin_visc  # Bubble Reynold's number

    dbl_p_col = math.tanh(math.sqrt(3 / 2 * (1 + (3 / 16 * dbl_re) / (1 + 0.249 * dbl_re ** 0.56))) * (particle_diam / dbl_bubble_diam)) ** 2  # Prob. collision, modified Luttrell and Yoon

    if dbl_p_col >= 1:
        dbl_p_col = 1

    dbl_eiw = GRAVITY / (4 * PI) * (WATER_VISCOSITY ** 3 / dbl_e_bulk) ** 0.25
    dbl_eka = (dbl_mass_bubble * dbl_u2_bulk - 2 * (dbl_bubble_diam / particle_diam) ** 2 * dbl_mass_particle * dbl_u1_bulk) ** 2 / (100 * (dbl_mass_bubble + 2 * (dbl_bubble_diam / particle_diam) ** 2 * dbl_mass_particle))

    dbl_p_i = 13 * math.sqrt((9 * WATER_VISCOSITY ** 2) / (dbl_bubble_diam * dbl_surface_tension * dbl_total_dens))
    dbl_pr = math.exp(-dbl_eiw / dbl_eka) if dbl_eka != 0 else 0
    dbl_pf_transfer = 1  # dbl_p_i * (1 - dbl_pr)

    # Froth Recovery - fitting parameter, can be optimized
    dbl_b = 2.3
    dbl_top_bubble_diam = dbl_bubble_diam * (1 / dbl_bbl_ratio)
    dbl_coverage_factor = 2 * (dbl_bubble_diam / particle_diam) ** 2
    dbl_buoyant_diam = dbl_bubble_diam * ((1 - 0.001275) / ((sg - 1) * dbl_coverage_factor)) ** (1 / 3)
    dbl_pfr = math.exp(dbl_b * particle_diam / dbl_buoyant_diam)
    dbl_rmax = dbl_bubble_diam / dbl_top_bubble_diam
    dbl_froth_ret_time = froth_height / sp_gas_rate * dbl_pfr
    dbl_alpha = 0.05

    dbl_r_attachment = dbl_rmax * math.exp(-dbl_alpha * dbl_froth_ret_time)

    dbl_qair = sp_gas_rate * dbl_cell_area
    dbl_qliq = cell_volume / (ret_time / num_cells * 60)
    dbl_r_water_max = (dbl_qair / dbl_qliq) / ((1 / 0.2) - 1)

    if dbl_r_water_max > 1:
        dbl_r_water_max = 0.1

    dbl_r_entrainment = dbl_r_water_max * math.exp(-0.0325 * (dbl_particle_dens - WATER_DENSITY) / 1000 - 0.063 * particle_diam * 1000000)

    dbl_froth_recovery_factor = dbl_r_entrainment + dbl_r_attachment

    # Cap froth recovery factor at 1.0 (cannot exceed 100%)
    if dbl_froth_recovery_factor > 1.0:
        dbl_froth_recovery_factor = 1.0

    # Rate Constant
    dbl_rate_const = dbl_beta * dbl_n_bubble * dbl_p_att * dbl_p_col * (1 - dbl_p_det) * 60  # x60 to make 1/min

    # Dispersion Model (uses Peclet number for axial dispersion)
    Aa = (1 + 4 * dbl_rate_const * ret_time / (num_cells * pe)) ** 0.5

    # Collection zone recovery using dispersion model (accounts for back-mixing)
    dbl_recovery_ci = 1 - 4 * Aa * math.exp(pe / 2) / ((1 + Aa) ** 2 * math.exp(Aa * pe / 2) - (1 - Aa) ** 2 * math.exp(-Aa * pe / 2))

    dbl_recovery_i = dbl_recovery_ci * dbl_froth_recovery_factor / (dbl_recovery_ci * dbl_froth_recovery_factor + 1 - dbl_recovery_ci)  # eq 6.2 finch & dobby

    if water_or_particle == "Water":
        return dbl_r_water_max
    elif water_or_particle == "Particle":
        return 1 - (1 - dbl_recovery_i) ** num_cells
    else:
        return 0.0


def col_flot_rec(sg, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
                 particle_z_pot, bubble_z_pot, num_cells, ret_time, cell_volume,
                 froth_height, frother, frother_conc, particle_diam, grade,
                 contact_angle, permitivity, dielectric, pe, water_or_particle,
                 dbl_cell_area, dbl_bbl_ratio):
    """
    Calculate column flotation recovery.
    Similar to flot_rec but with slightly different parameters for column flotation.
    """
    global dbl_energy_barrier, dbl_drag_beta

    # Convert particle diameter from microns to meters
    particle_diam = particle_diam * 0.000001

    # Constants
    dbl_h_c_factor = 150
    dbl_bulk_zone = 0.5
    dbl_impeller_zone = 15
    dbl_detach_f = 0.5  # Adjustable parameter for fitting
    dbl_particle_dens = sg * 1000  # x1000 for kg/m^3
    dbl_vol_imp_zone = 0.1  # Set impeller zone 1/10
    sp_power = sp_power * 1000  # x1000 for w/m^3
    sp_gas_rate = sp_gas_rate / 100  # /100 for m/s

    # Frother constants
    gamma_mibc = 0.000005  # mol/m^2
    gamma_ppg400 = 0.000001  # mol/m^2
    gamma_octanol = 0.000008  # mol/m^2
    gamma_pentanol = 0.000006  # mol/m^2
    k_mibc = 230  # M^-1
    k_ppg400 = 1700000  # M^-1
    k_octanol = 2200  # M^-1
    k_pentanol = 55  # M^-1

    # Froth parameters
    dbl_coverage = 0.3  # Different from flot_rec
    dbl_bubble_f = 0.5

    # Surface tension calculations
    if frother == 2:  # MIBC
        frother_conc = frother_conc / 102170
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_mibc * math.log(k_mibc * frother_conc + 1)
    elif frother == 3:  # PPG 400
        frother_conc = frother_conc / 134170
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_ppg400 * math.log(k_ppg400 * frother_conc + 1)
    elif frother == 4:  # Octanol
        frother_conc = frother_conc / 130230
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_octanol * math.log(k_octanol * frother_conc + 1)
    elif frother == 5:  # Pentanol
        frother_conc = frother_conc / 88150
        dbl_surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_pentanol * math.log(k_pentanol * frother_conc + 1)
    else:
        dbl_surface_tension = 0.07243

    # Energy Dissipation
    dbl_total_dens = air_fraction * AIR_DENSITY + (1 - air_fraction) * slurry_fraction * dbl_particle_dens + (1 - slurry_fraction) * WATER_DENSITY
    dbl_e_mean = sp_power / dbl_total_dens
    dbl_e_bulk = dbl_bulk_zone * dbl_e_mean
    dbl_e_impeller = dbl_impeller_zone * dbl_e_mean

    dbl_bubble_diam = dbl_bubble_f * (2.11 * dbl_surface_tension / (WATER_DENSITY * dbl_e_impeller ** 0.66)) ** 0.6

    dbl_num_attached = dbl_coverage * 4 * (dbl_bubble_diam / particle_diam) ** 2

    # Cell Calculations - note: collision diameter calculation is different from flot_rec
    dbl_collision_diam = particle_diam + dbl_bubble_diam  # Different from flot_rec (no 0.5 factor)
    dbl_vol_particle = (4 / 3) * PI * (particle_diam / 2) ** 3
    dbl_vol_bubble = (4 / 3) * PI * (dbl_bubble_diam / 2) ** 3
    dbl_vol_bp = dbl_vol_bubble + dbl_vol_particle
    dbl_kin_visc = WATER_VISCOSITY / WATER_DENSITY
    dbl_mass_particle = dbl_particle_dens * dbl_vol_particle
    dbl_mass_bubble = AIR_DENSITY * dbl_vol_bubble
    dbl_mass_bp = dbl_mass_bubble + dbl_mass_particle
    dbl_mass_total = cell_volume * dbl_total_dens

    # Velocities by Dissipation
    dbl_u1_bulk = (0.4 * (dbl_e_bulk ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (dbl_kin_visc ** (-1 / 3)) * (dbl_particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
    dbl_u2_bulk = 2 * (dbl_e_bulk * dbl_bubble_diam) ** (2 / 3)
    dbl_u1_mean = (0.4 * (dbl_e_mean ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (dbl_kin_visc ** (-1 / 3)) * (dbl_particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
    dbl_u2_mean = 2 * (dbl_e_mean * dbl_bubble_diam) ** (2 / 3)

    dbl_beta = (2 ** (3 / 2)) * (PI ** 0.5) * (dbl_collision_diam ** 2) * math.sqrt(dbl_u1_bulk + dbl_u2_bulk)

    # Calc # Density of Bubbles
    dbl_n_bubble = air_fraction / dbl_vol_bubble
    dbl_n_particle = (1 - air_fraction) * slurry_fraction / dbl_vol_particle
    dbl_z_bubb_particle = dbl_beta * dbl_n_bubble * dbl_n_particle

    dbl_work_adhesion = dbl_surface_tension * PI * (particle_diam / 2) ** 2 * (1 - math.cos(contact_angle * (PI / 180))) ** 2

    # Energy Barrier
    energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot)

    if dbl_energy_barrier <= 0:
        dbl_energy_barrier = 0

    # Kinetic Energy of Attachment
    dbl_kinetic_e_attach = 0.5 * dbl_mass_particle * dbl_u1_bulk / (dbl_drag_beta ** 2)
    dbl_kinetic_e_detach = 0.5 * dbl_mass_particle * (dbl_detach_f * (particle_diam + dbl_bubble_diam) * math.sqrt(dbl_e_impeller / dbl_kin_visc)) ** 2

    # Probabilities
    dbl_p_att = math.exp(-dbl_energy_barrier / dbl_kinetic_e_attach) if dbl_kinetic_e_attach != 0 else 0
    dbl_p_det = math.exp(-(dbl_work_adhesion + dbl_energy_barrier) / dbl_kinetic_e_detach) if dbl_kinetic_e_detach != 0 else 0
    dbl_re = math.sqrt(dbl_u2_bulk) * dbl_bubble_diam / dbl_kin_visc

    dbl_p_col = math.tanh(math.sqrt(3 / 2 * (1 + (3 / 16 * dbl_re) / (1 + 0.249 * dbl_re ** 0.56))) * (particle_diam / dbl_bubble_diam)) ** 2

    if dbl_p_col >= 1:
        dbl_p_col = 1

    dbl_eiw = GRAVITY / (4 * PI) * (WATER_VISCOSITY ** 3 / dbl_e_bulk) ** 0.25
    dbl_eka = (dbl_mass_bubble * dbl_u2_bulk - 2 * (dbl_bubble_diam / particle_diam) ** 2 * dbl_mass_particle * dbl_u1_bulk) ** 2 / (100 * (dbl_mass_bubble + 2 * (dbl_bubble_diam / particle_diam) ** 2 * dbl_mass_particle))

    dbl_p_i = 13 * math.sqrt((9 * WATER_VISCOSITY ** 2) / (dbl_bubble_diam * dbl_surface_tension * dbl_total_dens))
    dbl_pr = math.exp(-dbl_eiw / dbl_eka) if dbl_eka != 0 else 0
    dbl_pf_transfer = 1

    # Froth Recovery
    dbl_b = 2.2  # Different from flot_rec
    dbl_top_bubble_diam = dbl_bubble_diam * dbl_bbl_ratio  # Different from flot_rec (no 1/)
    dbl_coverage_factor = 2 * (dbl_bubble_diam / particle_diam) ** 2
    dbl_buoyant_diam = dbl_bubble_diam * ((1 - 0.001275) / ((sg - 1) * dbl_coverage_factor)) ** (1 / 3)
    dbl_pfr = math.exp(dbl_b * particle_diam / dbl_buoyant_diam)
    dbl_rmax = dbl_bubble_diam / dbl_top_bubble_diam
    dbl_froth_ret_time = froth_height / sp_gas_rate * dbl_pfr
    dbl_alpha = 0.05

    dbl_r_attachment = dbl_rmax * math.exp(-dbl_alpha * dbl_froth_ret_time)

    dbl_qair = sp_gas_rate * dbl_cell_area
    dbl_qliq = cell_volume / (ret_time / num_cells * 60)
    dbl_r_water_max = (dbl_qair / dbl_qliq) / ((1 / 0.2) - 1)

    if dbl_r_water_max > 1:
        dbl_r_water_max = 0.15  # Different from flot_rec

    dbl_r_entrainment = dbl_r_water_max * math.exp(-0.0325 * (dbl_particle_dens - WATER_DENSITY) / 1000 - 0.063 * particle_diam * 1000000)

    dbl_froth_recovery_factor = dbl_r_entrainment + dbl_r_attachment

    # Cap froth recovery factor at 1.0 (cannot exceed 100%)
    if dbl_froth_recovery_factor > 1.0:
        dbl_froth_recovery_factor = 1.0

    # Rate Constant
    dbl_rate_const = dbl_beta * dbl_n_bubble * dbl_p_att * dbl_p_col * (1 - dbl_p_det) * 60

    # Dispersion Model (uses Peclet number for axial dispersion)
    Aa = (1 + 4 * dbl_rate_const * ret_time / (num_cells * pe)) ** 0.5

    # Collection zone recovery using dispersion model (accounts for back-mixing)
    dbl_recovery_ci = 1 - 4 * Aa * math.exp(pe / 2) / ((1 + Aa) ** 2 * math.exp(Aa * pe / 2) - (1 - Aa) ** 2 * math.exp(-Aa * pe / 2))

    dbl_recovery_i = dbl_recovery_ci * dbl_froth_recovery_factor / (dbl_recovery_ci * dbl_froth_recovery_factor + 1 - dbl_recovery_ci)

    if water_or_particle == "Water":
        return dbl_r_water_max
    elif water_or_particle == "Particle":
        return 1 - (1 - dbl_recovery_i) ** num_cells
    else:
        return 0.0


def cyclone_rec(size, cut_size, alpha, water_rec):
    """
    Lynch-Rao Hydrocyclone model, with entrainment due to water recovery.

    Args:
        size: Particle size
        cut_size: Cut size (d50)
        alpha: Sharpness parameter
        water_rec: Water recovery

    Returns:
        float: Cyclone recovery
    """
    x = size / cut_size
    y = (math.exp(alpha * x) - 1) / (math.exp(alpha * x) + math.exp(alpha) - 2)
    cyclone_rec = y * (1 - water_rec) + water_rec

    return cyclone_rec


def grind_lib(tonnage_input, grade_input, lib_intensity):
    """
    Grinding model based on mass balance and liberation intensity.
    Liberation intensity indicates what portion of middling material is liberated.
    Liberation_intensity -> what degree of liberation you can achieve (fitting parameter).

    Args:
        tonnage_input: List/array of 4 tonnage values
        grade_input: List/array of 4 grade values
        lib_intensity: Liberation intensity (%)

    Returns:
        list: Updated tonnage values after liberation [4 elements]
    """
    tonnage = [0] + list(tonnage_input)  # Add dummy 0 at index 0 for 1-based indexing
    grade = [0] + list(grade_input)
    grind_lib_temp = [0] * 5  # Index 0-4, we'll use 1-4

    # Material in class = material in + material broken in - material broken out. Grade in = grade out
    grind_lib_temp[1] = tonnage[1] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[4]) / (grade[1] - grade[4]) + \
                        (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[4]) / (grade[1] - grade[4])
    grind_lib_temp[2] = tonnage[2] - tonnage[2] * lib_intensity / 100
    grind_lib_temp[3] = tonnage[3] - tonnage[3] * lib_intensity / 100
    grind_lib_temp[4] = tonnage[4] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[1]) / (grade[4] - grade[1]) + \
                        (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[1]) / (grade[4] - grade[1])

    return grind_lib_temp[1:5]  # Return elements 1-4


def grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity):
    """
    Grinding breakage model using population balance.

    Args:
        tonnage_input: List/array of 10 tonnage values
        top_size_input: List/array of 10 top size values
        bottom_size_input: List/array of 10 bottom size values
        selection_input: List/array of 10 selection function values
        break_intensity: Breakage intensity (number of cycles)

    Returns:
        numpy.ndarray: Updated tonnage distribution [10x1 array]
    """
    size_class = 10

    # Convert inputs to numpy arrays
    tonnage = np.zeros((11, 2))
    top_size = np.zeros((11, 2))
    bottom_size = np.zeros((11, 2))

    for i in range(1, 11):
        tonnage[i, 1] = tonnage_input[i-1]
        top_size[i, 1] = top_size_input[i-1]
        bottom_size[i, 1] = bottom_size_input[i-1]

        if i > 1:
            if tonnage[i, 1] == 0 and tonnage[i-1, 1] != 0:
                size_class = i - 1

    # Build breakage matrix
    breakage = np.zeros((11, 11))
    selection = np.zeros((11, 11))
    identity = np.zeros((11, 11))

    # Calculate breakage distribution
    for i in range(1, size_class + 1):
        if i == 1:
            breakage[i, 1] = 1 - (1 - math.exp(-bottom_size[i, 1] / top_size[1, 1])) / (1 - math.exp(-1))
        else:
            breakage[i, 1] = (1 - math.exp(-bottom_size[i-1, 1] / top_size[1, 1])) / (1 - math.exp(-1)) - \
                            (1 - math.exp(-bottom_size[i, 1] / top_size[1, 1])) / (1 - math.exp(-1))

    # Build full breakage and selection matrices
    for i in range(1, size_class + 1):
        column_sum = 0
        for j in range(1, size_class + 1):
            if i != 1 and j != size_class and j != 1:
                breakage[j, i] = breakage[j-1, i-1]

            if j == size_class:
                breakage[j, i] = 1 - column_sum
            else:
                column_sum = column_sum + breakage[j, i]

            if i == j:
                identity[i, j] = 1
                selection[i, j] = selection_input[i-1]
            else:
                identity[i, j] = 0
                selection[i, j] = 0

    # Matrix operations
    breakage_slice = breakage[1:size_class+1, 1:size_class+1]
    selection_slice = selection[1:size_class+1, 1:size_class+1]
    identity_slice = identity[1:size_class+1, 1:size_class+1]
    tonnage_slice = tonnage[1:size_class+1, 1:2]

    bs = np.matmul(breakage_slice, selection_slice)
    x = bs + identity_slice - selection_slice

    low_cycles = int(np.floor(break_intensity))
    high_cycles = low_cycles + 1

    # Calculate product for low cycles
    if low_cycles == 0:
        prod_low = tonnage_slice
    else:
        temp = x.copy()
        for i in range(low_cycles - 1):
            temp = np.matmul(temp, x)
        prod_low = np.matmul(temp, tonnage_slice)

    # Calculate product for high cycles
    if high_cycles == 0:
        prod_high = tonnage_slice
    else:
        temp = x.copy()
        for i in range(high_cycles - 1):
            temp = np.matmul(temp, x)
        prod_high = np.matmul(temp, tonnage_slice)

    # Interpolate between low and high cycles
    grind_break_temp = np.zeros((11, 2))
    for i in range(1, size_class + 1):
        grind_break_temp[i, 1] = prod_high[i-1, 0] * (break_intensity - low_cycles) + prod_low[i-1, 0] * (high_cycles - break_intensity)

    return grind_break_temp[1:11, 1]  # Return elements 1-10 from column 1

# tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity = 313, 75, 20, 

sg = 4.2 # Specific gravity
sp_power = 1 # Specific power (W/m^3 before conversion)
sp_gas_rate = 2 # Superficial gas rate (cm/s before conversion)
air_fraction = 0.05 # Air fraction
slurry_fraction = 0.15 # Slurry fraction
particle_z_pot = -50 # Particle zeta potential
bubble_z_pot = -.03 # Bubble zeta potential
num_cells = 10 # Number of cells in the flotation bank
ret_time = 23.67 # Total retention time (minutes) for entire bank
cell_volume = 700 #m^3
froth_height = .165 # Froth height (m)
frother = 4  # (2=MIBC, 3=PPG400, 4=Octanol, 5=Pentanol)
frother_conc = 50 # Frother concentration mg/L
particle_diam = 38.73 # Particle diameter (in microns)
grade = .2871 # Grade in % copper
contact_angle = 52.3 # Contact angle (degrees)
permitivity = 8.854*(10**-12) # Permitivity
dielectric = 86.5 # Dielectric constant
pe = 4 # Pe value 
water_or_particle = "Particle" # String "Water" or "Particle"
dbl_cell_area = 200 # Cell area
dbl_bbl_ratio = 10 # Bubble ratio

tonnage_input = 313 # Tonnage input
grade_input = 0.2871 # Grade input
lib_intensity = 50 # Liberation intensity (%)


# Size-by-grade flotation recovery calculation
# Calculate recovery for each size-grade combination

# Define size classes (microns) - from your matrix
size_classes = [212.13, 106.07, 38.73, 10.00]

# Define grade classes (% Chalco in film) and corresponding contact angles
grade_classes = [0.05, 0.20, 0.40, 0.75, 1.0]  # 0.1-10%, 10-30%, 30-50%, 50-100%, 100% Chalco
contact_angles = [16, 33, 42, 47, 50]  # degrees, corresponding to each grade class

# Mass distribution by size (total mass %)
total_mass_pct = [13.9, 15.8, 31.3, 39.2, 100.1]
micron_distribution = [212.13, 106.07, 38.73, 10.00]

throughput_grades = [0.05, 0.20, 0.40, 0.75, 1.0]  # 0.1-10%, 10-30%, 30-50%, 50-100%, 100% Chalco
throughput_distribution = [
    [136.22, .87, .42, .20, .29],
    [156.95, .61, .18, .14, .12],
    [309.24, .92, .61, .38, 1.85],
    [386.63, .76, .34, .39, 3.88],
    [989.04, 3.16, 1.55, 1.11, 6.14]
]

# Grade distribution within each size class (% of that size class)
grade_distribution = [
    [98.71, 0.63, 0.31, 0.14, 0.21],  # 212.13 micron
    [99.34, 0.38, 0.11, 0.09, 0.08],  # 106.07 micron
    [98.80, 0.30, 0.20, 0.12, 0.59],  # 38.73 micron
    [98.63, 0.19, 0.09, 0.10, 0.99],   # 10.00 micron
    [98.9, .32, .16, .11, .61] # average
]

sg_by_grade = [3.6, 4, 4.2, 4.2, 4.5]
contact_angles_by_grade = [10, 12.5, 27, 52.3, 60]
particle_z_pot_by_grade = [-50, -50, -50, -50, -50]
grade_by_grade = [0.002, 8, 12, 28, 34.65]


print("\n" + "="*100)
print("SIZE-BY-GRADE FLOTATION RECOVERY CALCULATION")
print("="*100)
print(f"{'Size (µm)':<12} {'Grade Range':<18} {'Mass (t)':<12} {'Grade %Cu':<12} {'Contact Angle':<15} {'Recovery %':<15}")
print("-"*100)

grade_labels = ["0.1-10%", "10-30%", "30-50%", "50-100%", "100%"]

total_feed_mass = 0
total_recovered_mass = 0

# Calculate recovery for each combination
for size_idx, particle_size in enumerate(size_classes):
    size_total_mass = total_mass_pct[size_idx]

    for grade_idx, grade_pct in enumerate(grade_classes):
        # Use grade-specific parameters
        contact_angle_deg = contact_angles_by_grade[grade_idx]
        sg_for_grade = sg_by_grade[grade_idx]
        particle_z_pot_for_grade = particle_z_pot_by_grade[grade_idx]
        grade_value = grade_by_grade[grade_idx]  # Actual grade % Cu

        grade_mass_fraction = grade_distribution[size_idx][grade_idx] / 100.0  # Convert % to fraction
        mass_in_class = throughput_distribution[size_idx][grade_idx]  # Use actual throughput in tonnes

        # Calculate recovery for this size-grade class
        recovery = flot_rec(
            sg_for_grade, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
            particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
            froth_height, frother, frother_conc, particle_size, grade_value,
            contact_angle_deg, permitivity, dielectric, pe, water_or_particle,
            dbl_cell_area, dbl_bbl_ratio
        )

        total_feed_mass += mass_in_class
        total_recovered_mass += mass_in_class * recovery

        print(f"{particle_size:<12.2f} {grade_labels[grade_idx]:<18} {mass_in_class:<12.2f} {grade_value:<12.2f} {contact_angle_deg:<15} {recovery*100:<15.2f}")
    print()  # Blank line between size classes

overall_recovery = total_recovered_mass / total_feed_mass if total_feed_mass > 0 else 0

print("="*100)
print(f"\nOVERALL WEIGHTED RECOVERY: {overall_recovery*100:.2f}%")
print(f"Total Feed Mass: {total_feed_mass:.2f} tonnes")
print(f"Total Recovered Mass: {total_recovered_mass:.2f} tonnes")
print(f"Total Cu in Feed: {sum([throughput_distribution[i][j] * grade_by_grade[j] / 100 for i in range(4) for j in range(5)]):.2f} tonnes Cu")
print(f"Total Cu Recovered: {sum([throughput_distribution[i][j] * grade_by_grade[j] / 100 * flot_rec(sg_by_grade[j], sp_power, sp_gas_rate, air_fraction, slurry_fraction, particle_z_pot_by_grade[j], bubble_z_pot, num_cells, ret_time, cell_volume, froth_height, frother, frother_conc, size_classes[i], grade_by_grade[j], contact_angles_by_grade[j], permitivity, dielectric, pe, water_or_particle, dbl_cell_area, dbl_bbl_ratio) for i in range(4) for j in range(5)]):.2f} tonnes Cu")
print("="*100)


tonnage_input = [138, 158, 313, 392, 100, 100, 100, 100, 100, 100]
top_size_input = [212.13, 106.07, 38.73, 10.00, 5.00, 2.50, 1.25, 0.625, 0.3125, 0.15625]
bottom_size_input = [106.07, 38.73, 10.00, 5.00, 2.50, 1.25, 0.625, 0.3125, 0.15625, 0.078125]
selection_input = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
break_intensity = 313
grinding_breaking = grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity)
print(grinding_breaking)