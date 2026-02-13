from PIL import Image, ImageDraw, ImageFont

# Paths to your images
image_files = [
    "graph5.png",  # D (top rectangle)
    "graph2.png",  # A (bottom left)
    "graph3.png",  # B (bottom middle)
    "graph4.png"   # C (bottom right)
]

# Load images
images = [Image.open(f) for f in image_files]

# ============================
# Set sizes for each image (width, height)
# ============================
size_top = (1300, 510)     # top rectangle
size_bottom = (500, 400)  # bottom squares

sizes = [size_top, size_bottom, size_bottom, size_bottom]

# Resize images according to sizes
images = [images[i].resize(sizes[i]) for i in range(4)]

# ============================
# Create canvas
# Width = max(top row width, bottom row width)
# Height = sum of top row + bottom row heights
# ============================
top_row_width = size_top[0]
bottom_row_width = size_bottom[0] * 3
canvas_width = max(top_row_width, bottom_row_width)
canvas_height = size_top[1] + size_bottom[1]

canvas = Image.new('RGB', (canvas_width, canvas_height), 'white')

# ============================
# Paste images
# ============================
# Top: D centered
top_position = ((canvas_width - size_top[0]) // 2, 0)
canvas.paste(images[0], top_position)

# Bottom row: B, C, D (centered horizontally)
total_bottom_width = size_bottom[0] * 3
start_x = (canvas_width - total_bottom_width) // 2
bottom_positions = [
    (start_x, size_top[1]),                  # bottom left
    (start_x + size_bottom[0], size_top[1]), # bottom middle
    (start_x + size_bottom[0]*2, size_top[1]) # bottom right
]
for i in range(1, 4):
    canvas.paste(images[i], bottom_positions[i-1])

# ============================
# Add labels (bold, bigger, in parentheses)
# ============================
draw = ImageDraw.Draw(canvas)

# Try to use a TrueType font for bigger bold letters
try:
    font = ImageFont.truetype("arialbd.ttf", 18)  # Arial Bold, size 48
except:
    font = ImageFont.load_default()  # fallback if arialbd.ttf is not available

labels = ["(a)", "(b)", "(c)", "(d)"]
# Top image label
draw.text((top_position[0]+10, top_position[1]+10), labels[0], fill="black", font=font)
# Bottom row labels
for i, pos in enumerate(bottom_positions):
    draw.text((pos[0]+10, pos[1]+10), labels[i+1], fill="black", font=font)

# ============================
# Save and show
# ============================
canvas.save("merged_collage.png")
canvas.show()
