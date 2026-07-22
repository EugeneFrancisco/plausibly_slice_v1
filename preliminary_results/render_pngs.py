"""Render the saved PLink projections as consistently sized PNG files."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from plink import LinkManager


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "png"
SIZE = 1600
MARGIN = 120
STROKE = 24
COLORS = ("#d62728", "#1f77b4", "#2ca02c")  # red, blue, green

DIAGRAMS = (
    ("friends/6_2_1_friend_K13n3596.lnk", "$6_2$ 1-friend: K13n3596"),
    (
        "friends/6_2_3_friend_K14n10164_mirror.lnk",
        "$6_2$ 3-friend: mirror of K14n10164",
    ),
    ("friends/conway_1_friend_A.lnk", "Conway 1-friend A"),
    ("friends/conway_1_friend_B.lnk", "Conway 1-friend B"),
    ("rbg_links/6_2_K13n3596_1_RBG.lnk", "$6_2$/K13n3596 1-RBG"),
    ("rbg_links/conway_friend_A_1_RBG.lnk", "Conway/friend A 1-RBG"),
)


def load_projection(path):
    manager = LinkManager()
    manager._from_string(path.read_text())
    return manager.polylines(break_at_overcrossings=True)


def render(path):
    components = load_projection(path)
    points = [point for lines, _ in components for line in lines for point in line]
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)
    scale = min(
        (SIZE - 2 * MARGIN) / max(max_x - min_x, 1),
        (SIZE - 2 * MARGIN) / max(max_y - min_y, 1),
    )
    width = (max_x - min_x) * scale
    height = (max_y - min_y) * scale
    offset_x = (SIZE - width) / 2 - min_x * scale
    offset_y = (SIZE - height) / 2 - min_y * scale

    def transform(point):
        return (point[0] * scale + offset_x, point[1] * scale + offset_y)

    image = Image.new("RGB", (SIZE, SIZE), "white")
    draw = ImageDraw.Draw(image)
    multi_component = len(components) == 3
    for index, (lines, _) in enumerate(components):
        color = COLORS[index] if multi_component else "#202020"
        for line in lines:
            transformed = [transform(point) for point in line]
            draw.line(
                transformed,
                fill=color,
                width=STROKE,
                joint="curve",
            )
            radius = STROKE / 2
            for x, y in (transformed[0], transformed[-1]):
                draw.ellipse(
                    (x - radius, y - radius, x + radius, y + radius),
                    fill=color,
                )
    return image


def make_contact_sheet(rendered):
    columns, rows = 3, 2
    cell_width, cell_height = 800, 820
    sheet = Image.new("RGB", (columns * cell_width, rows * cell_height), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default(size=32)
    for index, (image, label) in enumerate(rendered):
        thumbnail = image.resize((700, 700), Image.Resampling.LANCZOS)
        column, row = index % columns, index // columns
        x = column * cell_width + 50
        y = row * cell_height + 15
        sheet.paste(thumbnail, (x, y))
        clean_label = label.replace("$", "")
        box = draw.textbbox((0, 0), clean_label, font=font)
        text_width = box[2] - box[0]
        draw.text(
            (column * cell_width + (cell_width - text_width) / 2, y + 725),
            clean_label,
            fill="#202020",
            font=font,
        )
    return sheet


def main():
    OUTPUT.mkdir(exist_ok=True)
    rendered = []
    for relative_path, label in DIAGRAMS:
        source = HERE / relative_path
        image = render(source)
        destination = OUTPUT / f"{source.stem}.png"
        image.save(destination, dpi=(300, 300), optimize=True)
        rendered.append((image, label))
        print(destination.name, image.size)
    contact_sheet = make_contact_sheet(rendered)
    contact_sheet.save(OUTPUT / "all_six_diagrams.png", dpi=(300, 300), optimize=True)
    print("all_six_diagrams.png", contact_sheet.size)


if __name__ == "__main__":
    main()
