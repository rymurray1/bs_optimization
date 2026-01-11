import math
import numpy as np

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


# Test grind_break
if __name__ == "__main__":
    # Test with multiple size classes (more realistic)
    # Size classes from coarse to fine
    tonnage_input = [100, 150, 50, 10, 3]  # tonnes in each class
    top_size_input = [600, 300, 150, 75, 37.5]  # microns (top of each class)
    bottom_size_input = [300, 150, 75, 37.5, 18.75]  # microns (bottom of each class)
    selection_input = [0.8, 0.6, 0.4, 0.2, 0.1]  # coarse particles break more easily
    break_intensity = 5

    print("GRINDING TEST - MULTI-SIZE CLASS")
    print("=" * 80)
    print(f"Break intensity: {break_intensity}")
    print()
    print("INPUT FEED DISTRIBUTION:")
    print(f"{'Size Class':>12} | {'Size Range (um)':>20} | {'Tonnage':>10} | {'Selection':>10}")
    print("-" * 80)

    total_feed = sum(tonnage_input)
    for i, (tonnage, top, bottom, sel) in enumerate(zip(tonnage_input, top_size_input, bottom_size_input, selection_input)):
        pct = tonnage / total_feed * 100
        print(f"{i+1:12d} | {bottom:8.2f} - {top:8.2f} | {tonnage:10.1f} | {sel:10.2f}")

    print("-" * 80)
    print(f"{'TOTAL':>12} | {'':>20} | {total_feed:10.1f} |")
    print()

    # Run grinding model
    ground_tonnage = grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity)

    print("OUTPUT PRODUCT DISTRIBUTION:")
    print(f"{'Size Class':>12} | {'Size Range (um)':>20} | {'Tonnage':>10} | {'% of Total':>12}")
    print("-" * 80)

    total_product = np.sum(ground_tonnage)
    for i, (tonnage, top, bottom) in enumerate(zip(ground_tonnage, top_size_input, bottom_size_input)):
        pct = tonnage / total_product * 100
        print(f"{i+1:12d} | {bottom:8.2f} - {top:8.2f} | {tonnage:10.2f} | {pct:11.2f}%")

    print("-" * 80)
    print(f"{'TOTAL':>12} | {'':>20} | {total_product:10.2f} |")
    print()

    print("SUMMARY:")
    print(f"Total feed: {total_feed:.2f} tonnes")
    print(f"Total product: {total_product:.2f} tonnes")
    print(f"Mass balance check: {abs(total_feed - total_product) < 0.01}")
    print()

    # Show how material redistributed
    print("MASS REDISTRIBUTION:")
    print(f"{'Size Class':>12} | {'Feed':>10} | {'Product':>10} | {'Change':>10}")
    print("-" * 80)
    for i, (feed, product) in enumerate(zip(tonnage_input, ground_tonnage)):
        change = product - feed
        print(f"{i+1:12d} | {feed:10.1f} | {product:10.2f} | {change:+10.2f}")
    print()

    print("INTERPRETATION:")
    print("- Positive change: Material generated from breakage of coarser particles")
    print("- Negative change: Material broken into finer size classes")
    print("- Coarse classes typically show negative change (material breaking down)")
    print("- Fine classes typically show positive change (receiving broken material)")