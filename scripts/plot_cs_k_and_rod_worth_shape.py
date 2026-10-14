from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cs_prediction_validation_slurm_5615" / "cs_k_and_rod_worth_shape.png"

Z = [64.0, 72.0, 80.0, 83.638895, 88.0, 92.5]
K = [1.086778179, 1.059766925, 1.021845315, 1.000016390, 0.969979640, 0.937435656]
M = [-0.003013356296, -0.0040968633712, -0.0055917107824,
     -0.0062914608992, -0.0068244248264, -0.0061670356952]


def font(size, bold=False):
    name = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path("C:/Windows/Fonts") / name), size)


def hermite(az, ak, am, bz, bk, bm, t):
    h = bz - az
    t2, t3 = t * t, t * t * t
    value = ((2*t3-3*t2+1)*ak + (t3-2*t2+t)*h*am
             + (-2*t3+3*t2)*bk + (t3-t2)*h*bm)
    derivative = (((6*t2-6*t)/h)*ak + (3*t2-4*t+1)*am
                  + ((-6*t2+6*t)/h)*bk + (3*t2-2*t)*bm)
    return value, derivative


curve = []
for i in range(len(Z)-1):
    steps = 240
    stop = steps + 1 if i == len(Z)-2 else steps
    for s in range(stop):
        t = s / steps
        kval, slope = hermite(Z[i], K[i], M[i], Z[i+1], K[i+1], M[i+1], t)
        zval = Z[i] + t * (Z[i+1] - Z[i])
        worth = -slope / (kval*kval) * 1e5
        curve.append((zval, kval, slope, worth))

critical_z = None
for a, b in zip(curve, curve[1:]):
    if (a[1]-1)*(b[1]-1) <= 0:
        q = (1-a[1])/(b[1]-a[1])
        critical_z = a[0] + q*(b[0]-a[0])
        break
peak = max(curve, key=lambda p: p[3])

W, H = 1800, 1320
img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)
blue, orange, grid, text, gray = "#1769aa", "#d55e00", "#d7dce2", "#17202a", "#59636e"

d.text((90, 38), "CS point values and local slopes: reconstructed rod-worth shape", fill=text, font=font(42, True))
d.line((92, 103, 140, 103), fill=blue, width=6)
d.text((153, 87), "Cubic Hermite curve", fill=text, font=font(23))
d.ellipse((445, 94, 461, 110), fill=blue)
d.text((475, 87), "Calculated CS point", fill=text, font=font(23))
d.line((760, 103, 808, 103), fill=orange, width=5)
d.text((821, 87), "Local CS slope", fill=text, font=font(23))


def panel(box, ymin, ymax, ylabel, title, yticks):
    x0, y0, x1, y1 = box
    def sx(x): return x0 + (x-Z[0])/(Z[-1]-Z[0])*(x1-x0)
    def sy(y): return y1 - (y-ymin)/(ymax-ymin)*(y1-y0)
    d.text((x0, y0-48), title, fill=text, font=font(28, True))
    for val in yticks:
        yy = sy(val)
        d.line((x0, yy, x1, yy), fill=grid, width=2)
        label = f"{val:.3f}" if ymax < 2 else f"{val:.0f}"
        bbox = d.textbbox((0,0), label, font=font(20))
        d.text((x0-16-(bbox[2]-bbox[0]), yy-(bbox[3]-bbox[1])/2), label, fill=gray, font=font(20))
    for val in [64, 68, 72, 76, 80, 84, 88, 92]:
        xx = sx(val)
        d.line((xx, y0, xx, y1), fill=grid, width=2)
        label = str(val)
        bbox = d.textbbox((0,0), label, font=font(20))
        d.text((xx-(bbox[2]-bbox[0])/2, y1+12), label, fill=gray, font=font(20))
    d.rectangle(box, outline="#89939e", width=2)
    d.text(((x0+x1)//2-145, y1+48), "Rod insertion depth, z (cm)", fill=text, font=font(23))
    d.text((25, (y0+y1)//2-14), ylabel, fill=text, font=font(23))
    return sx, sy


top = (150, 190, 1735, 670)
sx, sy = panel(top, 0.925, 1.100, "k_eff", "Reconstructed k_eff(z)", [0.94, 0.98, 1.02, 1.06, 1.10])
d.line([(sx(p[0]), sy(p[1])) for p in curve], fill=blue, width=6, joint="curve")
for zi, ki, mi in zip(Z, K, M):
    lo, hi = max(Z[0], zi-1.25), min(Z[-1], zi+1.25)
    d.line((sx(lo), sy(ki+mi*(lo-zi)), sx(hi), sy(ki+mi*(hi-zi))), fill=orange, width=5)
    x, y = sx(zi), sy(ki)
    d.ellipse((x-9, y-9, x+9, y+9), fill=blue, outline="white", width=3)
d.line((top[0], sy(1.0), top[2], sy(1.0)), fill=gray, width=3)
if critical_z is not None:
    cx = sx(critical_z)
    d.line((cx, top[1], cx, top[3]), fill="#7d8792", width=3)
    d.text((cx+14, sy(1.0)-44), f"k = 1 at about {critical_z:.3f} cm", fill=text, font=font(22))

bottom = (150, 805, 1735, 1190)
sx2, sy2 = panel(bottom, 150, 900, "-drho/dz (pcm/cm)", "Differential rod worth implied by the reconstructed curve", [200, 350, 500, 650, 800])
poly = [(sx2(curve[0][0]), bottom[3])] + [(sx2(p[0]), sy2(p[3])) for p in curve] + [(sx2(curve[-1][0]), bottom[3])]
d.polygon(poly, fill="#e6f0f8")
d.line([(sx2(p[0]), sy2(p[3])) for p in curve], fill=blue, width=6, joint="curve")
for zi, ki, mi in zip(Z, K, M):
    worth = -mi/(ki*ki)*1e5
    x, y = sx2(zi), sy2(worth)
    d.ellipse((x-9, y-9, x+9, y+9), fill=blue, outline="white", width=3)
px, py = sx2(peak[0]), sy2(peak[3])
d.ellipse((px-10, py-10, px+10, py+10), fill=orange)
d.line((px-5, py+14, px-120, py+82), fill=gray, width=3)
d.text((px-410, py+76), f"maximum about {peak[3]:.0f} pcm/cm at {peak[0]:.2f} cm", fill=text, font=font(22))

d.text((150, 1268), "Interpolation only between calculated depths. The 92.5 cm endpoint slope uses L10 and may retain source-tail bias.", fill=gray, font=font(20))
OUT.parent.mkdir(parents=True, exist_ok=True)
img.save(OUT)
print(OUT)
