#!/usr/bin/env python3
#
# make_radiator_view_json.py - interactively mark the lower-left and
# upper-right corners of the diamond crystal in a radiator view PNG,
# enter its physical width/height in mm, and write out a sidecar JSON
# file describing the crystal's geometry within the image, for use by
# spotfinder.py's get_radiator_geometry().
#
# This is a standalone developer tool, run once per new radiator image
# on a local workstation -- it is not part of the spotfinder web app
# itself, and its dependency on matplotlib does not affect spotfinder.py.
#
# author: richard.t.jones at uconn.edu
#
# usage: python3 make_radiator_view_json.py JD80-211_front_view.png

import sys
import os
import json
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <radiator_view.png>")
        sys.exit(1)

    png_path = sys.argv[1]
    if not os.path.isfile(png_path):
        print(f"error: file not found: {png_path}")
        sys.exit(1)

    img = mpimg.imread(png_path)
    height_px_full, width_px_full = img.shape[0], img.shape[1]

    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.set_title("Click the LOWER-LEFT corner of the radiator,\n"
                 "then the UPPER-RIGHT corner")
    plt.tight_layout()

    print(f"Image native size: {width_px_full} x {height_px_full} px")
    print("Click the LOWER-LEFT corner of the radiator crystal, "
          "then the UPPER-RIGHT corner, in the image window...")
    pts = plt.ginput(2, timeout=0)
    plt.close(fig)

    if len(pts) != 2:
        print("error: did not receive two clicks, aborting")
        sys.exit(1)

    (x0, y0), (x1, y1) = pts
    # matplotlib imshow reports (x, y) in native image pixel coordinates,
    # with y=0 at the top row -- same orientation as the raw pixel grid.
    # min/max are used here rather than trusting click order, so the
    # result is robust even if the two corners are clicked in a different
    # order than "lower-left then upper-right".
    x_left = min(x0, x1)
    x_right = max(x0, x1)
    y_top = min(y0, y1)
    y_bottom = max(y0, y1)

    width_px = x_right - x_left
    height_px = y_bottom - y_top
    center_x_px = (x_left + x_right) / 2
    center_y_px = (y_top + y_bottom) / 2

    print(f"Selected region: {width_px:.1f} x {height_px:.1f} px, "
          f"center at ({center_x_px:.1f}, {center_y_px:.1f})")

    width_mm = float(input("Enter the physical width of the radiator (mm): "))
    height_mm = float(input("Enter the physical height of the radiator (mm): "))

    geom = {
        "width_mm": width_mm,
        "height_mm": height_mm,
        "width_px": round(width_px),
        "height_px": round(height_px),
        "center_x_px": round(center_x_px),
        "center_y_px": round(center_y_px),
    }

    json_path = os.path.splitext(png_path)[0] + ".json"
    with open(json_path, "w") as f:
        json.dump(geom, f, indent=2)

    print(f"Wrote {json_path}:")
    print(json.dumps(geom, indent=2))


if __name__ == "__main__":
    main()
