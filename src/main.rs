extern crate sdl2;
extern crate time;

// implementation of http://www.karlsims.com/rd.html

use sdl2::pixels::PixelFormatEnum;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;

use time::PreciseTime;

const HEIGHT: usize = 640;
const WIDTH: usize = 640;
const GRID_HEIGHT: usize = 200;
const GRID_WIDTH: usize = 200;

type FloatingPoint = f64;

const INIT_DA: f64 = 1.0;
const INIT_DB: f64 = 0.5;
// const INIT_F: f64 = 0.0367;
// const INIT_K: f64 = 0.0649;
const INIT_F: f64 = 0.055;
const INIT_K: f64 = 0.062;
const INIT_DELTAT: f64 = 1.0;
const INIT_SIM_SPEED: f64 = 1.0;

const LAPLACIAN_IDX_OFFSETS: [(isize, isize); 9] =
    [(-1, -1), (0, -1), (1, -1),
     (-1, 0), (0, 0), (1, 0),
     (-1, 1), (0, 1), (1, 1)];

const W_CENTER: f64 = -1.0;
const W_ADJACENT: f64 = 0.2;
const W_DIAGONAL: f64 = 0.05;

const LAPLACIAN_WEIGHTS: [FloatingPoint; 9] =
    [W_DIAGONAL, W_ADJACENT, W_DIAGONAL,
     W_ADJACENT,   W_CENTER, W_ADJACENT,
     W_DIAGONAL, W_ADJACENT, W_DIAGONAL];

const GRID_A: usize = 0;
const GRID_B: usize = 1;

struct ReactionDiffusionGrid {
    width: usize,
    height: usize,
    d_a: FloatingPoint,
    d_b: FloatingPoint,
    f: FloatingPoint,
    k: FloatingPoint,
    grids: [Box<[FloatingPoint]>; 2],
}

impl ReactionDiffusionGrid {
    fn new(width: usize, height: usize, dA: FloatingPoint, dB: FloatingPoint, f: FloatingPoint, k: FloatingPoint) -> Self {
        let mut vecA = vec![0.0; width * height];
        let mut vecB = vec![0.0; width * height];

        return ReactionDiffusionGrid { width: width, height: height, d_a: dA, d_b: dB, f: f, k: k,
            grids: [vecA.into_boxed_slice(), vecB.into_boxed_slice()] }
    }

    fn seed(&mut self, i:usize, x: usize, y: usize, size_x: usize, size_y: usize, value: FloatingPoint) {
        for off_y in 0..size_y {
            for off_x in 0..size_x {
                let relem = self.grid_val(i, (x + off_x) as isize, (y + off_y) as isize);
                *relem = value;
            }
        }
    }


    fn seed_round(&mut self, i:usize, x: usize, y: usize, size: usize, value: FloatingPoint) {
        let isz = size as isize;
        for off_y in -isz..isz {
            for off_x in -isz..isz {
                if off_x * off_x + off_y * off_y > isz * isz {
                    continue;
                }
                let relem = self.grid_val(i, x as isize + off_x, y as isize + off_y);
                *relem = value;
            }
        }
    }

    fn update(&mut self, deltaT: FloatingPoint) {
        let mut tmp_grid = vec![0.0; self.width * self.height];

        for i in 0..2 {
            for y in 0..(self.height as isize) {
                for x in 0..(self.width as isize) {
                    let idx = self.idx(x, y);
                    let a = self.grids[GRID_A][idx];
                    let b = self.grids[GRID_B][idx];
                    let abb = a * b * b;
                    let laplacian = self.laplacian(i, x, y);
                    // if laplacian == 0 {
                    //     continue;
                    // }
                    if i == GRID_A {
                        // println!("{:?}", laplacian);

                        let mut val =
                            self.grids[GRID_A][idx] + (self.d_a * laplacian - abb + self.f * (1.0 - a)) * deltaT;
                        // if val > 1.0 || val < 0.0 {
                        //     panic!("Illegal val for A! {} {} {} {}", val, self.grids[GRID_A][idx], a, b);
                        // }
                        val = softclamp(val);
                        tmp_grid[idx] = val;
                    } else {
                        let mut val =
                            self.grids[GRID_B][idx] + (self.d_b * laplacian + abb - (self.k + self.f) * b) * deltaT;
                        // if val > 1.0 || val < 0.0 {
                        //     panic!("Illegal val for B! {} {} {} {}", val, self.grids[GRID_B][idx], a, b);
                        // }
                        val = softclamp(val);
                        tmp_grid[idx] = val;
                    }
                }
            }

            self.grids[i].copy_from_slice(tmp_grid.as_slice());
        }
    }

    fn laplacian(&mut self, i: usize, x: isize, y: isize) -> FloatingPoint {
        let mut off_x: isize = 0;
        let mut off_y: isize = 0;
        let mut weight: FloatingPoint = 0.0;
        let mut laplacian = 0.0;
        for k in 0..9 {
            let offsets = LAPLACIAN_IDX_OFFSETS[k];
            off_x = offsets.0;
            off_y = offsets.1;
            weight = LAPLACIAN_WEIGHTS[k];
            let ngb_val = *self.grid_val(i, x + off_x, y + off_y);

            laplacian += weight * ngb_val;
        }
        return laplacian;
    }

    fn idx(&self, x: isize, y: isize) -> usize {
        let mut nx = x;
        let mut ny = y;
        let iwidth = self.width as isize;
        let iheight = self.height as isize;
        if nx < 0 { nx += iwidth; }
        else if nx >= iwidth { nx -= iwidth; }
        if ny < 0 { ny += iheight; }
        else if ny >= iheight { ny -= iheight; }
        return (ny * iwidth + nx) as usize
        // (modulus(y, self.height) * self.width + modulus(x, self.width)) as usize
    }

    fn grid_val(&mut self, i: usize, x: isize, y: isize) -> &mut FloatingPoint {
        return &mut self.grids[i][self.idx(x, y)];
    }

    fn debug_dump(&self) {
        // println!("{:?}", self.grids[0]);
        println!("{:?}", self.grids[1]);
    }
}

fn modulus(a: isize, b: usize) -> usize {
    let ib = b as isize;
    (((a % ib) + ib) % ib) as usize
}

fn softclamp(f: FloatingPoint) -> FloatingPoint {
    if f < -100.0 { return -100.0; }
    else if f > 100.0 { return 100.0; }
    else { return f; }
}

pub fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("rust-sdl2 demo: Video",
        HEIGHT as u32, WIDTH as u32)
    .position_centered()
    .opengl()
    .build()
    .unwrap();

    let mut renderer = window.renderer().build().unwrap();

    let mut texture = renderer.create_texture_streaming(
        PixelFormatEnum::RGB24, GRID_HEIGHT as u32, GRID_WIDTH as u32).unwrap();

    let mut event_pump = sdl_context.event_pump().unwrap();


    let mut grid = ReactionDiffusionGrid::new(GRID_WIDTH, GRID_HEIGHT, INIT_DA, INIT_DB, INIT_F, INIT_K);

    grid.seed(GRID_A, 0, 0, GRID_WIDTH, GRID_HEIGHT, 1.0);
    grid.seed_round(GRID_B, GRID_WIDTH / 2 - 35, GRID_HEIGHT / 2 - 4, GRID_WIDTH / 50, 1.0);
    grid.seed_round(GRID_B, GRID_WIDTH / 2 + 15, GRID_HEIGHT / 2 + 30, GRID_WIDTH / 50, 0.7);

    let start = PreciseTime::now();

    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..}
                | Event::KeyDown { keycode: Some(Keycode::Escape), .. }
                | Event::KeyDown { keycode: Some(Keycode::Q), .. } => {
                    break 'running
                }
                ,
                Event::KeyDown { keycode: Some(Keycode::D), .. } => {
                    grid.debug_dump();
                },
                _ => {}
            }
        }


        let t = PreciseTime::now();
        texture.with_lock(None, |buffer: &mut [u8], pitch: usize| {
            for y in 0..GRID_HEIGHT {
                for x in 0..GRID_WIDTH {
                    let offset = y * pitch + x * 3;
                    let val_A = *grid.grid_val(GRID_A, x as isize, y as isize);
                    let val_B = *grid.grid_val(GRID_B, x as isize, y as isize);

                    // let g = (0.5 * val_A + 0.5 * val_B) * 255.0;
                    // buffer[offset + 0] = g as u8;
                    // buffer[offset + 1] = g as u8;
                    // buffer[offset + 2] = g as u8;
                    buffer[offset + 0] = 0*(val_B * 255.0) as u8;
                    buffer[offset + 1] = 1*(val_B * 255.0) as u8;
                    buffer[offset + 2] = 0*(val_B * 255.0) as u8;
                }
            }
        }).unwrap();

        renderer.clear();
        renderer.copy(&texture, None, None);
        renderer.present();
        // println!("Frame time: {}", t.to(PreciseTime::now()).num_milliseconds());

        grid.update(INIT_DELTAT);
    }
}
