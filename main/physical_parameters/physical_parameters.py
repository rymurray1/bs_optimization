import math
import numpy as np

# Constants (these would typically be defined elsewhere or passed as parameters)
PI = math.pi
WATER_DENSITY = 1000  # kg/m^3
AIR_DENSITY = 1.225  # kg/m^3
WATER_VISCOSITY = 0.001  # Pa·s
GRAVITY = 9.81  # m/s^2

# Global variables for energy barrier calculation
energy_barrier_value = 0.0
drag_beta_value = 0.0


def calc_energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot):
    """
    Calculate energy barrier for particle-bubble attachment.
    """
    global energy_barrier_value, drag_beta_value
    # TODO: Implement proper energy barrier calculation from VBA
    # TEMPORARY: Using calibrated value to achieve ~90% recovery target
    # This should be replaced with actual DLVO calculation
    # For copper flotation, energy barriers should be very low for valuable minerals
    energy_barrier_value = 1e-18  # Calibrated value in Joules (targeting ~90% overall recovery)
    drag_beta_value = 1.0


def flot_rec(physical_variables, operating_conditions, cell_design, model_parameters, physical_constants, control_parameters):
    """
    Calculate flotation recovery.

    Args:
        specific_gravity: Specific gravity
        specific_power_input: Specific power (W/m^3 before conversion)
        specific_gas_rate: Specific gas rate (cm/s before conversion)
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
        cell_area: Cell area
        bbl_ratio: Bubble ratio

    Returns:
        float: Flotation recovery value
    """
    global energy_barrier_value, drag_beta_value
    
    #physical variables
    specific_gravity, particle_diam, particle_z_pot, bubble_z_pot, contact_angle, grade = physical_variables
    
    #operating conditions
    specific_power_input, specific_gas_velocity, air_fraction, slurry_fraction, frother_conc, froth_type = operating_conditions
    
    #cell design
    num_cells, ret_time, cell_volume, froth_height, cell_area, bbl_ratio = cell_design
    
    #model parameters
    b, alpha, coverage, bubble_f, detach_f, bulk_zone, pe = model_parameters
    
    #physical constants
    permitivity, dielectric = physical_constants
    
    #control parameters
    water_or_particle = control_parameters

    b = b
    alpha = alpha
    coverage = coverage # Froth parameters
    bubble_f = bubble_f # Froth parameters
    detach_f = detach_f  # Adjustable parameter for fitting - varies depending on ore; hydrodynamic parameters; optimized using ML
    bulk_zone = bulk_zone # Can be changed, varies depending on ore and hydrodynamic parameters    
    # Convert particle diameter from microns to meters
    particle_diam = particle_diam * 0.000001

    # Constants
    h_c_factor = 150
    impeller_zone = 15
    particle_dens = specific_gravity * 1000  # x1000 for kg/m^3
    vol_imp_zone = 0.1  # Set impeller zone 1/10
    specific_power_input = specific_power_input* 1000  # x1000 for w/m^3
    specific_gas_rate = specific_gas_velocity / 100  # /100 for m/s

    # Frother constants
    gamma_mibc = 0.000005  # mol/m^2
    gamma_ppg400 = 0.000001  # mol/m^2
    gamma_octanol = 0.000008  # mol/m^2
    gamma_pentanol = 0.000006  # mol/m^2
    k_mibc = 230  # M^-1
    k_ppg400 = 1700000  # M^-1
    k_octanol = 2200  # M^-1
    k_pentanol = 55  # M^-1

    # Surface tension calculations
    if froth_type== 2:  # MIBC
        frother_conc = frother_conc / 102170  # Convert ppm to mol/L
        surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_mibc * math.log(k_mibc * frother_conc + 1)
    elif froth_type== 3:  # PPG 400
        frother_conc = frother_conc / 134170  # Convert ppm to mol/L
        surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_ppg400 * math.log(k_ppg400 * frother_conc + 1)
    elif froth_type== 4:  # Octanol
        frother_conc = frother_conc / 130230  # Convert ppm to mol/L
        surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_octanol * math.log(k_octanol * frother_conc + 1)
    elif froth_type== 5:  # Pentanol
        frother_conc = frother_conc / 88150  # Convert ppm to mol/L
        surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_pentanol * math.log(k_pentanol * frother_conc + 1)
    else:
        surface_tension = 0.07243  # Pure water @ 23°C

    # Energy Dissipation
    total_dens = air_fraction * AIR_DENSITY + (1 - air_fraction) * slurry_fraction * particle_dens + (1 - slurry_fraction) * WATER_DENSITY
    e_mean = specific_power_input/ total_dens
    e_bulk = bulk_zone * e_mean
    e_impeller = impeller_zone * e_mean

    bubble_diam = bubble_f * (2.11 * surface_tension / (WATER_DENSITY * e_impeller ** 0.66)) ** 0.6

    num_attached = coverage * 4 * (bubble_diam / particle_diam) ** 2  # Num of particles attached to one bubble

    # Cell Calculations
    collision_diam = 0.5 * (particle_diam + bubble_diam)  # Avg diam of collision
    vol_particle = (4 / 3) * PI * (particle_diam / 2) ** 3  # Vol 1 part.
    vol_bubble = (4 / 3) * PI * (bubble_diam / 2) ** 3  # Vol 1 bubb.
    vol_bp = vol_bubble + vol_particle  # Vol of 1 BP aggregate
    kin_visc = WATER_VISCOSITY / WATER_DENSITY
    mass_particle = particle_dens * vol_particle  # Mass 1 part.
    mass_bubble = AIR_DENSITY * vol_bubble  # Mass 1 bubb.
    mass_bp = mass_bubble + mass_particle  # Mass of 1 BP aggregate
    mass_total = cell_volume * total_dens

    # Velocities by Dissipation
    u1_bulk = (0.4 * (e_bulk ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (kin_visc ** (-1 / 3)) * (particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2  # For attachment
    u2_bulk = 2 * (e_bulk * bubble_diam) ** (2 / 3)
    u1_mean = (0.4 * (e_mean ** (4 / 9)) * (particle_diam ** (7 / 9)) *
                   (kin_visc ** (-1 / 3)) * (particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
    u2_mean = 2 * (e_mean * bubble_diam) ** (2 / 3)

    beta = (2 ** (3 / 2)) * (PI ** 0.5) * (collision_diam ** 2) * math.sqrt(u1_bulk + u2_bulk)  # From Abrahamson model using bulk dissipation

    # Calc # Density of Bubbles
    n_bubble = air_fraction / vol_bubble
    n_particle = (1 - air_fraction) * slurry_fraction / vol_particle
    z_bubb_particle = beta * n_bubble * n_particle

    work_adhesion = surface_tension * PI * (particle_diam / 2) ** 2 * (1 - math.cos(contact_angle * (PI / 180))) ** 2  # Calc work of adhesion for 1 particle

    # Energy Barrier
    calc_energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot)

    if energy_barrier_value <= 0:
        energy_barrier_value = 0

    # Kinetic Energy of Attachment
    kinetic_e_attach = 0.5 * mass_particle * u1_bulk / (drag_beta_value ** 2)
    kinetic_e_detach = 0.5 * mass_particle * (detach_f * (particle_diam + bubble_diam) * math.sqrt(e_impeller / kin_visc)) ** 2

    # Probabilities
    p_att = math.exp(-energy_barrier_value / kinetic_e_attach) if kinetic_e_attach != 0 else 0  # Prob. of attachment
    p_det = math.exp(-(work_adhesion + energy_barrier_value) / kinetic_e_detach) if kinetic_e_detach != 0 else 0  # Prob. of detachment
    re = math.sqrt(u2_bulk) * bubble_diam / kin_visc  # Bubble Reynold's number

    p_col = math.tanh(math.sqrt(3 / 2 * (1 + (3 / 16 * re) / (1 + 0.249 * re ** 0.56))) * (particle_diam / bubble_diam)) ** 2  # Prob. collision, modified Luttrell and Yoon

    if p_col >= 1:
        p_col = 1

    eiw = GRAVITY / (4 * PI) * (WATER_VISCOSITY ** 3 / e_bulk) ** 0.25
    eka = (mass_bubble * u2_bulk - 2 * (bubble_diam / particle_diam) ** 2 * mass_particle * u1_bulk) ** 2 / (100 * (mass_bubble + 2 * (bubble_diam / particle_diam) ** 2 * mass_particle))

    p_i = 13 * math.sqrt((9 * WATER_VISCOSITY ** 2) / (bubble_diam * surface_tension * total_dens))
    pr = math.exp(-eiw / eka) if eka != 0 else 0
    pf_transfer = 1  # p_i * (1 - pr)

    top_bubble_diam = bubble_diam * (1 / bbl_ratio)
    coverage_factor = 2 * (bubble_diam / particle_diam) ** 2
    buoyant_diam = bubble_diam * ((1 - 0.001275) / ((specific_gravity - 1) * coverage_factor)) ** (1 / 3)
    pfr = math.exp(b * particle_diam / buoyant_diam)
    rmax = bubble_diam / top_bubble_diam
    froth_ret_time = froth_height / specific_gas_rate * pfr
    
    r_attachment = rmax * math.exp(-alpha * froth_ret_time)

    air_flow_rate = specific_gas_rate * cell_area
    water_flow_rate = cell_volume / (ret_time / num_cells * 60)
    r_water_max = (air_flow_rate / water_flow_rate) / ((1 / 0.2) - 1)

    if r_water_max > 1:
        r_water_max = 0.1

    r_entrainment = r_water_max * math.exp(-0.0325 * (particle_dens - WATER_DENSITY) / 1000 - 0.063 * particle_diam * 1000000)

    froth_recovery_factor = r_entrainment + r_attachment


    if froth_recovery_factor > 1.0:
        froth_recovery_factor = 1.0

    # Rate Constant
    rate_const = beta * n_bubble * p_att * p_col * (1 - p_det) * 60  # x60 to make 1/min

    # Dispersion Model (uses Peclet number for axial dispersion)
    Aa = (1 + 4 * rate_const * ret_time / (num_cells * pe)) ** 0.5

    # Collection zone recovery using dispersion model (accounts for back-mixing)
    recovery_ci = 1 - 4 * Aa * math.exp(pe / 2) / ((1 + Aa) ** 2 * math.exp(Aa * pe / 2) - (1 - Aa) ** 2 * math.exp(-Aa * pe / 2))

    recovery_i = recovery_ci * froth_recovery_factor / (recovery_ci * froth_recovery_factor + 1 - recovery_ci)  # eq 6.2 finch & dobby

    if water_or_particle == "Water":
        return r_water_max
    elif water_or_particle == "Particle":
        return 1 - (1 - recovery_i) ** num_cells
    else:
        return 0.0

# physical_variables = specific_gravity, particle_diam, particle_z_pot, bubble_z_pot, contact_angle, grade
# operating_conditions = specific_power_input, specific_gas_velocity, air_fraction, slurry_fraction, frother_conc, froth_type
# cell_design = num_cells, ret_time, cell_volume, froth_height, cell_area, bbl_ratio 
# model_parameters = b, alpha, coverage, bubble_f, detach_f, bulk_zone, pe
# physical_constants = permitivity, dielectric
# control_parameters = water_or_particle


# def col_flot_rec(specific_gravity, specific_power_input, specific_gas_rate, air_fraction, slurry_fraction,
#                  particle_z_pot, bubble_z_pot, num_cells, ret_time, cell_volume,
#                  froth_height, frother, frother_conc, particle_diam, grade,
#                  contact_angle, permitivity, dielectric, pe, water_or_particle,
#                  cell_area, bbl_ratio):
#     """
#     Calculate column flotation recovery.
#     Similar to flot_rec but with slightly different parameters for column flotation.
#     """
#     global energy_barrier_value, drag_beta_value

#     # Convert particle diameter from microns to meters
#     particle_diam = particle_diam * 0.000001

#     # Constants
#     h_c_factor = 150
#     detach_f = .55
#     bulk_zone = 0.5
#     impeller_zone = 15
#     particle_dens = specific_gravity * 1000  # x1000 for kg/m^3
#     vol_imp_zone = 0.1  # Set impeller zone 1/10
#     specific_power_input = specific_power_input* 1000  # x1000 for w/m^3
#     specific_gas_rate = specific_gas_rate / 100  # /100 for m/s

#     # Frother constants
#     gamma_mibc = 0.000005  # mol/m^2
#     gamma_ppg400 = 0.000001  # mol/m^2
#     gamma_octanol = 0.000008  # mol/m^2
#     gamma_pentanol = 0.000006  # mol/m^2
#     k_mibc = 230  # M^-1
#     k_ppg400 = 1700000  # M^-1
#     k_octanol = 2200  # M^-1
#     k_pentanol = 55  # M^-1

#     # Froth parameters
#     coverage = 0.3  # Different from flot_rec
#     bubble_f = 0.5

#     # Surface tension calculations
#     if froth_type== 2:  # MIBC
#         frother_conc = frother_conc / 102170
#         surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_mibc * math.log(k_mibc * frother_conc + 1)
#     elif froth_type== 3:  # PPG 400
#         frother_conc = frother_conc / 134170
#         surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_ppg400 * math.log(k_ppg400 * frother_conc + 1)
#     elif froth_type== 4:  # Octanol
#         frother_conc = frother_conc / 130230
#         surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_octanol * math.log(k_octanol * frother_conc + 1)
#     elif froth_type== 5:  # Pentanol
#         frother_conc = frother_conc / 88150
#         surface_tension = 0.07243 - 8.314 * (273.15 + 23) * gamma_pentanol * math.log(k_pentanol * frother_conc + 1)
#     else:
#         surface_tension = 0.07243

#     # Energy Dissipation
#     total_dens = air_fraction * AIR_DENSITY + (1 - air_fraction) * slurry_fraction * particle_dens + (1 - slurry_fraction) * WATER_DENSITY
#     e_mean = specific_power_input/ total_dens
#     e_bulk = bulk_zone * e_mean
#     e_impeller = impeller_zone * e_mean

#     bubble_diam = bubble_f * (2.11 * surface_tension / (WATER_DENSITY * e_impeller ** 0.66)) ** 0.6

#     num_attached = coverage * 4 * (bubble_diam / particle_diam) ** 2

#     # Cell Calculations - note: collision diameter calculation is different from flot_rec
#     collision_diam = particle_diam + bubble_diam  # Different from flot_rec (no 0.5 factor)
#     vol_particle = (4 / 3) * PI * (particle_diam / 2) ** 3
#     vol_bubble = (4 / 3) * PI * (bubble_diam / 2) ** 3
#     vol_bp = vol_bubble + vol_particle
#     kin_visc = WATER_VISCOSITY / WATER_DENSITY
#     mass_particle = particle_dens * vol_particle
#     mass_bubble = AIR_DENSITY * vol_bubble
#     mass_bp = mass_bubble + mass_particle
#     mass_total = cell_volume * total_dens

#     # Velocities by Dissipation
#     u1_bulk = (0.4 * (e_bulk ** (4 / 9)) * (particle_diam ** (7 / 9)) *
#                    (kin_visc ** (-1 / 3)) * (particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
#     u2_bulk = 2 * (e_bulk * bubble_diam) ** (2 / 3)
#     u1_mean = (0.4 * (e_mean ** (4 / 9)) * (particle_diam ** (7 / 9)) *
#                    (kin_visc ** (-1 / 3)) * (particle_dens / WATER_DENSITY - 1) ** (2 / 3)) ** 2
#     u2_mean = 2 * (e_mean * bubble_diam) ** (2 / 3)

#     beta = (2 ** (3 / 2)) * (PI ** 0.5) * (collision_diam ** 2) * math.sqrt(u1_bulk + u2_bulk)

#     # Calc # Density of Bubbles
#     n_bubble = air_fraction / vol_bubble
#     n_particle = (1 - air_fraction) * slurry_fraction / vol_particle
#     z_bubb_particle = beta * n_bubble * n_particle

#     work_adhesion = surface_tension * PI * (particle_diam / 2) ** 2 * (1 - math.cos(contact_angle * (PI / 180))) ** 2

#     # Energy Barrier
#     calc_energy_barrier(particle_diam, dielectric, permitivity, contact_angle, bubble_z_pot, particle_z_pot)

#     if energy_barrier_value <= 0:
#         energy_barrier_value = 0

#     # Kinetic Energy of Attachment
#     kinetic_e_attach = 0.5 * mass_particle * u1_bulk / (drag_beta_value ** 2)
#     kinetic_e_detach = 0.5 * mass_particle * (detach_f * (particle_diam + bubble_diam) * math.sqrt(e_impeller / kin_visc)) ** 2

#     # Probabilities
#     p_att = math.exp(-energy_barrier_value / kinetic_e_attach) if kinetic_e_attach != 0 else 0
#     p_det = math.exp(-(work_adhesion + energy_barrier_value) / kinetic_e_detach) if kinetic_e_detach != 0 else 0
#     re = math.sqrt(u2_bulk) * bubble_diam / kin_visc

#     p_col = math.tanh(math.sqrt(3 / 2 * (1 + (3 / 16 * re) / (1 + 0.249 * re ** 0.56))) * (particle_diam / bubble_diam)) ** 2

#     if p_col >= 1:
#         p_col = 1

#     eiw = GRAVITY / (4 * PI) * (WATER_VISCOSITY ** 3 / e_bulk) ** 0.25
#     eka = (mass_bubble * u2_bulk - 2 * (bubble_diam / particle_diam) ** 2 * mass_particle * u1_bulk) ** 2 / (100 * (mass_bubble + 2 * (bubble_diam / particle_diam) ** 2 * mass_particle))

#     p_i = 13 * math.sqrt((9 * WATER_VISCOSITY ** 2) / (bubble_diam * surface_tension * total_dens))
#     pr = math.exp(-eiw / eka) if eka != 0 else 0
#     pf_transfer = 1

#     # Froth Recovery
#     b = 2.2  # Different from flot_rec
#     top_bubble_diam = bubble_diam * bbl_ratio  # Different from flot_rec (no 1/)
#     coverage_factor = 2 * (bubble_diam / particle_diam) ** 2
#     buoyant_diam = bubble_diam * ((1 - 0.001275) / ((specific_gravity - 1) * coverage_factor)) ** (1 / 3)
#     pfr = math.exp(b * particle_diam / buoyant_diam)
#     rmax = bubble_diam / top_bubble_diam
#     froth_ret_time = froth_height / specific_gas_rate * pfr
#     alpha = 0.05

#     r_attachment = rmax * math.exp(-alpha * froth_ret_time)

#     air_flow_rate = specific_gas_rate * cell_area
#     water_flow_rate = cell_volume / (ret_time / num_cells * 60)
#     r_water_max = (air_flow_rate / water_flow_rate) / ((1 / 0.2) - 1)

#     if r_water_max > 1:
#         r_water_max = 0.15  # Different from flot_rec

#     r_entrainment = r_water_max * math.exp(-0.0325 * (particle_dens - WATER_DENSITY) / 1000 - 0.063 * particle_diam * 1000000)

#     froth_recovery_factor = r_entrainment + r_attachment

#     # Cap froth recovery factor at 1.0 (cannot exceed 100%)
#     if froth_recovery_factor > 1.0:
#         froth_recovery_factor = 1.0

#     # Rate Constant
#     rate_const = beta * n_bubble * p_att * p_col * (1 - p_det) * 60

#     # Dispersion Model (uses Peclet number for axial dispersion)
#     Aa = (1 + 4 * rate_const * ret_time / (num_cells * pe)) ** 0.5

#     # Collection zone recovery using dispersion model (accounts for back-mixing)
#     recovery_ci = 1 - 4 * Aa * math.exp(pe / 2) / ((1 + Aa) ** 2 * math.exp(Aa * pe / 2) - (1 - Aa) ** 2 * math.exp(-Aa * pe / 2))

#     recovery_i = recovery_ci * froth_recovery_factor / (recovery_ci * froth_recovery_factor + 1 - recovery_ci)

#     if water_or_particle == "Water":
#         return r_water_max
#     elif water_or_particle == "Particle":
#         return 1 - (1 - recovery_i) ** num_cells
#     else:
#         return 0.0

# def cyclone_rec(size, cut_size, alpha, water_rec):
#     """
#     Lynch-Rao Hydrocyclone model, with entrainment due to water recovery.

#     Args:
#         size: Particle size
#         cut_size: Cut size (d50)
#         alpha: Sharpness parameter
#         water_rec: Water recovery

#     Returns:
#         float: Cyclone recovery
#     """
#     x = size / cut_size
#     y = (math.exp(alpha * x) - 1) / (math.exp(alpha * x) + math.exp(alpha) - 2)
#     cyclone_rec = y * (1 - water_rec) + water_rec

#     return cyclone_rec


# # def grind_lib(tonnage_input, grade_input, lib_intensity):
#     """
#     Grinding model based on mass balance and liberation intensity.
#     Liberation intensity indicates what portion of middling material is liberated.
#     Liberation_intensity -> what degree of liberation you can achieve (fitting parameter).

#     Args:
#         tonnage_input: List/array of 4 tonnage values
#         grade_input: List/array of 4 grade values
#         lib_intensity: Liberation intensity (%)

#     Returns:
#         list: Updated tonnage values after liberation [4 elements]
#     """
#     tonnage = [0] + list(tonnage_input)  # Add dummy 0 at index 0 for 1-based indexing
#     grade = [0] + list(grade_input)
#     grind_lib_temp = [0] * 5  # Index 0-4, we'll use 1-4

#     # Material in class = material in + material broken in - material broken out. Grade in = grade out
#     grind_lib_temp[1] = tonnage[1] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[4]) / (grade[1] - grade[4]) + \
#                         (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[4]) / (grade[1] - grade[4])
#     grind_lib_temp[2] = tonnage[2] - tonnage[2] * lib_intensity / 100
#     grind_lib_temp[3] = tonnage[3] - tonnage[3] * lib_intensity / 100
#     grind_lib_temp[4] = tonnage[4] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[1]) / (grade[4] - grade[1]) + \
#                         (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[1]) / (grade[4] - grade[1])

#     return grind_lib_temp[1:5]  # Return elements 1-4

# def grind_lib(tonnage_input, grade_input, lib_intensity):
#     """
#     Grinding model based on mass balance and liberation intensity.
#     Liberation intensity indicates what portion of middling material is liberated.
#     Liberation_intensity -> what degree of liberation you can achieve (fitting parameter).

#     Args:
#         tonnage_input: List/array of 4 tonnage values
#         grade_input: List/array of 4 grade values
#         lib_intensity: Liberation intensity (%)

#     Returns:
#         list: Updated tonnage values after liberation [4 elements]
#     """
#     tonnage = [0] + list(tonnage_input)  # Add dummy 0 at index 0 for 1-based indexing
#     grade = [0] + list(grade_input)
#     grind_lib_temp = [0] * 5  # Index 0-4, we'll use 1-4

#     # Material in class = material in + material broken in - material broken out. Grade in = grade out
#     grind_lib_temp[1] = tonnage[1] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[4]) / (grade[1] - grade[4]) + \
#                         (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[4]) / (grade[1] - grade[4])
#     grind_lib_temp[2] = tonnage[2] - tonnage[2] * lib_intensity / 100
#     grind_lib_temp[3] = tonnage[3] - tonnage[3] * lib_intensity / 100
#     grind_lib_temp[4] = tonnage[4] + (lib_intensity / 100 * tonnage[2] * grade[2] - lib_intensity / 100 * tonnage[2] * grade[1]) / (grade[4] - grade[1]) + \
#                         (lib_intensity / 100 * tonnage[3] * grade[3] - lib_intensity / 100 * tonnage[3] * grade[1]) / (grade[4] - grade[1])

#     return grind_lib_temp[1:5]  # Return elements 1-4

def grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity):
    """
    Grinding breakage model using population balance.

    Args:
        tonnage_input: List/array of tonnage values (flexible length)
        top_size_input: List/array of top size values (same length as tonnage_input)
        bottom_size_input: List/array of bottom size values (same length as tonnage_input)
        selection_input: List/array of selection function values (same length as tonnage_input)
        break_intensity: Breakage intensity (number of cycles)

    Returns:
        numpy.ndarray: Updated tonnage distribution
    """
    # Determine size_class from input length
    num_size_classes = len(tonnage_input)
    size_class = num_size_classes

    # Convert inputs to numpy arrays (pad to allow 1-based indexing)
    tonnage = np.zeros((num_size_classes + 1, 2))
    top_size = np.zeros((num_size_classes + 1, 2))
    bottom_size = np.zeros((num_size_classes + 1, 2))

    for i in range(1, num_size_classes + 1):
        tonnage[i, 1] = tonnage_input[i-1]
        top_size[i, 1] = top_size_input[i-1]
        bottom_size[i, 1] = bottom_size_input[i-1]

        if i > 1:
            if tonnage[i, 1] == 0 and tonnage[i-1, 1] != 0:
                size_class = i - 1

    # Build breakage matrix
    breakage = np.zeros((num_size_classes + 1, num_size_classes + 1))
    selection = np.zeros((num_size_classes + 1, num_size_classes + 1))
    identity = np.zeros((num_size_classes + 1, num_size_classes + 1))

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
    grind_break_temp = np.zeros((num_size_classes + 1, 2))
    for i in range(1, size_class + 1):
        grind_break_temp[i, 1] = prod_high[i-1, 0] * (break_intensity - low_cycles) + prod_low[i-1, 0] * (high_cycles - break_intensity)

    return grind_break_temp[1:num_size_classes+1, 1]  # Return all size classes from column 1

def copper_leaching_rate(Ea, T, P_O2, n, Phi, H_plus, MFeS2, X_Cu, k0=1.0, R=8.314):
    """
    Calculate the copper leaching rate based on shrinking core kinetics.

    Equation: dX_Cu/dt = 3*k0*exp(-Ea/RT)*(P_O2_eff)^n*(1-X_Cu)^(2/3)

    Args:
        Ea: Activation energy (kJ/mol)
        T: Temperature (°C)
        P_O2: Partial pressure of oxygen (bar)
        n: Stoichiometry parameter (dimensionless)
        Phi: Mass of solids per solution volume (kg/L)
        H_plus: H+ ion concentration, typically expressed as -log(pH)
        MFeS2: Mass fraction of FeS2 (dimensionless)
        X_Cu: Molar ratio of Cu converted (fraction, 0-1)
        k0: Pre-exponential factor (1/s), default=1.0
        R: Universal gas constant (J/mol·K), default=8.314

    Returns:
        float: Rate of copper conversion (dX_Cu/dt) in 1/s

    Base case parameters (from image):
        Ea: 70-90 kJ/mol
        T: 180-202°C
        P_O2: 10-15 bar
        n: 0.5-1
        Phi: 10 kg/L
        H_plus: -log(pH)
        MFeS2: varies
    """
    # Convert temperature from Celsius to Kelvin
    T_K = T + 273.15

    # Convert Ea from kJ/mol to J/mol for consistency with R
    Ea_J = Ea * 1000

    # Calculate effective oxygen pressure (P_O2_eff = P_O2 / (1 + alpha*MFeS2))
    # For now, assuming P_O2_eff = P_O2 based on the equation shown
    # This can be modified if the relationship P_O2/(1+alpha*MFeS2) needs to be applied
    P_O2_eff = P_O2

    # Calculate the leaching rate
    # dX_Cu/dt = 3*k0*exp(-Ea/RT)*(P_O2_eff)^n*(1-X_Cu)^(2/3)
    exponential_term = math.exp(-Ea_J / (R * T_K))
    pressure_term = P_O2_eff ** n
    conversion_term = (1 - X_Cu) ** (2/3)

    dX_Cu_dt = 3 * k0 * exponential_term * pressure_term * conversion_term

    return dX_Cu_dt

def solvent_extraction(C_Cu_in_aq, K_ex, RH, H_plus, O_A):
    """
    Calculate output copper concentration from solvent extraction.

    Equation: C_Cu_out^aq = C_Cu_in^aq / (1 + K_ex * [RH]^2 / [H+]^2 * (O/A))

    Args:
        C_Cu_in_aq: Amount of Cu in feed solution to SX extractant (moles of Cu)
        K_ex: Extraction equilibrium constant (dimensionless)
        RH: Extractant concentration (molar)
        H_plus: Hydrogen ion concentration (molar)
        O_A: Organic/aqueous solution ratio (dimensionless)

    Returns:
        float: Output copper concentration (C_Cu_out^aq) in moles
    """
    # Calculate the denominator term
    denominator = 1 + K_ex * (RH ** 2) / (H_plus ** 2) * O_A

    # Calculate output concentration
    C_Cu_out_aq = C_Cu_in_aq / denominator

    return C_Cu_out_aq

def electrowinning_concentration(C_Cu_in, eta_CE, I, n, F, M_Cu):
    """
    Calculate output copper concentration from electrowinning.

    Equation: C_Cu_out = C_Cu_in - (eta_CE * I) / (n * F * F)

    Args:
        C_Cu_in: Input copper concentration (g/mol or appropriate units)
        eta_CE: Current efficiency (dimensionless, 0-1)
        I: Current (Amperes)
        n: Number of electrons per Cu (2 e- per Cu)
        F: Faraday constant (96485 C/mol)
        M_Cu: Molecular weight of Cu (63.55 g/mol)

    Returns:
        float: Output copper concentration (C_Cu_out)
    """
    # Calculate concentration change due to electrowinning
    C_Cu_out = C_Cu_in - (eta_CE * I) / (n * F * F)

    return C_Cu_out

def cell_voltage(E_eq, eta_cath, eta_anode, I, R):
    """
    Calculate cell voltage for electrowinning.

    Equation: V_cell = E_eq + eta_cath + eta_anode + IR

    Args:
        E_eq: Nernst equilibrium potential (V)
        eta_cath: Overpotentials at cathode (V)
        eta_anode: Overpotentials at anode (V)
        I: Current (A)
        R: Resistance (Ohms)

    Returns:
        float: Cell voltage (V)
    """
    V_cell = E_eq + eta_cath + eta_anode + I * R

    return V_cell

def nernst_potential(E_Cu_0, a_Cu2, R=8.314, T=298.15, F=96485, n=2):
    """
    Calculate Nernst equilibrium potential.

    Equation: E_eq = E_Cu^0 + (RT/2F) * ln(a_Cu2+)

    Args:
        E_Cu_0: Standard potential for Cu (V), typically 0.34 V
        a_Cu2: Activity of Cu2+ ions (dimensionless)
        R: Universal gas constant (J/mol·K), default=8.314
        T: Temperature (K), default=298.15
        F: Faraday constant (C/mol), default=96485
        n: Number of electrons (2 for Cu), default=2

    Returns:
        float: Nernst equilibrium potential (V)
    """
    E_eq = E_Cu_0 + (R * T / (n * F)) * math.log(a_Cu2)

    return E_eq

def overpotential(a, j, b):
    """
    Calculate overpotential as function of current density.

    Equation: eta = a * ln(j) + b

    Args:
        a: Tafel slope coefficient (V)
        j: Current density (I/A)
        b: Constant (V)

    Returns:
        float: Overpotential (V)
    """
    eta = a * math.log(j) + b

    return eta

# Base case parameters for electrowinning
ELECTROWINNING_BASE_PARAMS = {
    'M_Cu': 63.55,  # g/mol
    'n': 2,  # electrons per Cu
    'F': 96485,  # C/mol
    'eta_CE': 0.95,  # Current efficiency (midpoint of 0.9-0.98)
    'E_Cu_0': 0.34,  # Standard potential for Cu (V)
    'eta_cath_range': (0.1, 0.3),  # Overpotentials cathode (V)
    'eta_anode_range': (0.3, 0.6),  # Overpotentials anode (V)
    'T': 298.15,  # Temperature (K)
    'R': 8.314  # Universal gas constant (J/mol·K)
}

# Base case parameters for copper leaching
LEACHING_BASE_PARAMS = {
    'Ea': 80.0,  # Activation energy (kJ/mol) - midpoint of 70-90 range
    'T': 191.0,  # Temperature (°C) - midpoint of 180-202 range
    'P_O2': 12.5,  # Partial oxygen pressure (bar) - midpoint of 10-15 range
    'n': 0.75,  # Stoichiometry - midpoint of 0.5-1 range
    'Phi': 10.0,  # Mass of solids/solution volume (kg/L)
    'H_plus': 1.0,  # -log(pH), pH = 1 typical for acidic leaching
    'MFeS2': 0.1,  # Mass fraction of FeS2 - example value
    'X_Cu_initial': 0.0,  # Initial copper conversion (0 = no conversion yet)
    'k0': 1.0,  # Pre-exponential factor (1/s)
    'R': 8.314  # Universal gas constant (J/mol·K)
}

def flot_constants():
    return {"specific_gravity": 4.2,
            "specific_power_input": 1,
            "specific_gas_rate": 2,
            "air_fraction": 0.05,
            "slurry_fraction": 0.15,
            "particle_z_pot": -50,
            "bubble_z_pot": -0.03,
            "num_cells": 1,
            "ret_time": 23.67,
            "cell_volume": 700,
            "froth_height": 0.165,
            "froth_type": 4,
            "frother_conc": 50,
            "particle_diam": 75,
            "grade": 28,
            "contact_angle": 52.3,
            "permitivity": 8.854,
            "dielectric": 86.5,
            "pe": 4,
            "water_or_particle": "Particle",
            "cell_area": 200,
            "bbl_ratio": 10
            }
specific_gravity = 4.2 # Specific gravity
specific_power_input= 1 # Specific power (W/m^3 before conversion)
specific_gas_rate = 2 # specified gas rate (cm/s before conversion)
air_fraction = 0.05 # Air fraction
slurry_fraction = 0.15 # Slurry fraction
particle_z_pot = -50 # Particle zeta potential
bubble_z_pot = -.03 # Bubble zeta potential
num_cells = 1 # Number of cells in the flotation bank
ret_time = 23.67 # Total retention time (minutes) for entire bank
cell_volume = 700 #m^3
froth_height = .165 # Froth height (m)
froth_type = 4  # (2=MIBC, 3=PPG400, 4=Octanol, 5=Pentanol)
frother_conc = 50 # Frother concentration mg/L
particle_diam = 75 # Particle diameter (in microns)
grade = 28 # Grade in % copper
contact_angle = 52.3 # Contact angle (degrees)
permitivity = 8.854 # Permitivity
dielectric = 86.5 # Dielectric constant
pe = 4 # Pe value 
water_or_particle = "Particle" # String "Water" or "Particle"
cell_area = 200 # Cell area
bbl_ratio = 10 # Bubble ratio


def fitting_parameters():
    return {
        'b': 2,
        'alpha': 0.10,
        'coverage': 0.525,
        'bubble_f': 0.825,
        'detach_f': 0.5,
        'bulk_zone': 0.5
    }

def input_variables(): 
    return {
        "tonnage_input": 313,
        "grade_input": 28,
        "lib_intensity": 100,
        "specific_gravity_for_grade": 4.2,
        "particle_z_pot_for_grade": -50,
        "grade_value": 29,
        "particle_size": 75,
        "contact_angle_deg": 52.3
    }

b, alpha, coverage, bubble_f, detach_f, bulk_zone = fitting_parameters().values() 
tonnage_input, grade_input, lib_intensity, specific_gravity_for_grade, particle_z_pot_for_grade, grade_value, particle_size, contact_angle_deg = input_variables().values()  

recovery = flot_rec(
            specific_gravity_for_grade, specific_power_input, specific_gas_rate, air_fraction, slurry_fraction,
            particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
            froth_height, froth_type, frother_conc, particle_size, grade_value,
            contact_angle_deg, permitivity, dielectric, pe, water_or_particle,
            cell_area, bbl_ratio, b, alpha, coverage, bubble_f, detach_f, bulk_zone
        )
print(recovery)
# Size-by-grade flotation recovery calculation
# Calculate recovery for each size-grade combination

# # Define size classes (microns) - from your matrix
# size_classes = [150, 75]

# # Define grade classes (% Chalco in film) and corresponding contact angles
# grade_classes = [0.20, 0.40, 0.75]  # 0.1-10%, 10-30%, 30-50%, 50-100%, 100% Chalco
# contact_angles = [33, 42, 47]  # degrees, corresponding to each grade class

# # Mass distribution by size (total mass %)
# total_mass_pct = [15.8, 31.3, 39.2]
# micron_distribution = size_classes

# throughput_grades = grade_classes   # 0.1-10%, 10-30%, 30-50%, 50-100%, 100% Chalco
# throughput_distribution = [
#     [136.22, .87, .42, .20, .29],
#     [156.95, .61, .18, .14, .12],
#     [309.24, .92, .61, .38, 1.85],
#     [386.63, .76, .34, .39, 3.88],
#     [989.04, 3.16, 1.55, 1.11, 6.14]
# ]

# # Grade distribution within each size class (% of that size class)
# grade_distribution = [
#     [98.71, 0.63, 0.31, 0.14, 0.21],  # 212.13 micron
#     [99.34, 0.38, 0.11, 0.09, 0.08],  # 106.07 micron
#     [98.80, 0.30, 0.20, 0.12, 0.59],  # 38.73 micron
#     [98.63, 0.19, 0.09, 0.10, 0.99],   # 10.00 micron
#     [98.9, .32, .16, .11, .61] # average
# ]

# specific_gravity_by_grade = [3.6, 4, 4.2, 4.2, 4.5]
# contact_angles_by_grade = [22, 32.5, 42, 51.5, 62.5]  # Tuned to achieve ~3% concentrate grade
# particle_z_pot_by_grade = [-50, -50, -50, -50, -50]
# grade_by_grade = [0.002, 8, 12, 28, 34.65]


# print("\n" + "="*100)
# print("SIZE-BY-GRADE FLOTATION RECOVERY CALCULATION")
# print("="*100)
# print(f"{'Size (µm)':<12} {'Grade Range':<18} {'Mass (t)':<12} {'Grade %Cu':<12} {'Contact Angle':<15} {'Recovery %':<15}")
# print("-"*100)

# grade_labels = ["10-30%", "30-50%", "50-100%"]

# total_feed_mass = 0
# total_recovered_mass = 0
# total_cu_in_feed = 0
# total_cu_recovered = 0

# # Store concentrate data for each class
# concentrate_data = []

# # Calculate recovery for each combination
# for size_idx, particle_size in enumerate(size_classes):
#     size_total_mass = total_mass_pct[size_idx]

#     for grade_idx, grade_pct in enumerate(grade_classes):
#         # Use grade-specific parameters
#         contact_angle_deg = contact_angles_by_grade[grade_idx]
#         specific_gravity_for_grade = specific_gravity_by_grade[grade_idx]
#         particle_z_pot_for_grade = particle_z_pot_by_grade[grade_idx]
#         grade_value = grade_by_grade[grade_idx]  # Actual grade % Cu

#         grade_mass_fraction = grade_distribution[size_idx][grade_idx] / 100.0  # Convert % to fraction
#         mass_in_class = throughput_distribution[size_idx][grade_idx]  # Use actual throughput in tonnes

#         # Calculate recovery for this size-grade class
#         recovery = flot_rec(
#             specific_gravity_for_grade, specific_power_input, specific_gas_rate, air_fraction, slurry_fraction,
#             particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
#             froth_height, froth_type, frother_conc, particle_size, grade_value,
#             contact_angle_deg, permitivity, dielectric, pe, water_or_particle,
#             cell_area, bbl_ratio, b, alpha, coverage, bubble_f, detach_f, bulk_zone
#         )

#         # Calculate Cu content
#         cu_in_feed = mass_in_class * grade_value / 100
#         cu_recovered = cu_in_feed * recovery
#         mass_recovered = mass_in_class * recovery

#         total_feed_mass += mass_in_class
#         total_recovered_mass += mass_recovered
#         total_cu_in_feed += cu_in_feed
#         total_cu_recovered += cu_recovered

#         # Store concentrate data
#         concentrate_data.append({
#             'size': particle_size,
#             'size_idx': size_idx,
#             'grade_label': grade_labels[grade_idx],
#             'grade_idx': grade_idx,
#             'mass_feed': mass_in_class,
#             'mass_conc': mass_recovered,
#             'cu_feed': cu_in_feed,
#             'cu_conc': cu_recovered,
#             'grade_value': grade_value,
#             'recovery': recovery
#         })

#         print(f"{particle_size:<12.2f} {grade_labels[grade_idx]:<18} {mass_in_class:<12.2f} {grade_valuevalue:<12.2f} {contact_angle_deg:<15} {recovery*100:<15.2f}")
#     print()  # Blank line between size classes

# overall_recovery = total_recovered_mass / total_feed_mass if total_feed_mass > 0 else 0
# overall_concentrate_grade = (total_cu_recovered / total_recovered_mass * 100) if total_recovered_mass > 0 else 0

# print("="*100)
# print(f"\nOVERALL WEIGHTED RECOVERY: {overall_recovery*100:.2f}%")
# print(f"Total Feed Mass: {total_feed_mass:.2f} tonnes")
# print(f"Total Recovered Mass: {total_recovered_mass:.2f} tonnes")
# print(f"Total Cu in Feed: {total_cu_in_feed:.2f} tonnes Cu")
# print(f"Total Cu Recovered: {total_cu_recovered:.2f} tonnes Cu")
# print(f"Overall Concentrate Grade: {overall_concentrate_grade:.2f}% Cu")
# print("="*100)

tonnage_input = [1000,10000]
top_size_input = [150, 75]
bottom_size_input = [75, 20]
selection_input = [0.2, 0.3]
break_intensity = 8
grinding_breaking = grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity)
print(grinding_breaking)

"""
All fitting parameters:
bulk_zone = 0.5
impeller_zone = 15
detach_f = 0.55
coverage = 0.5
bubble_f = 0.725
b = 2.3
alpha = 0.1
energy_barrier_value (line 25) = 1e-18
"""
