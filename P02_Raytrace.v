module main
import os
import math
import gfx

////////////////////////////////////////////////////////////////////////////////////////
// Comment out lines in array below to prevent re-rendering every scene.
// If you create a new scene file, add it to the list below.
//
// NOTE: **BEFORE** you submit your solution, uncomment all lines, so
//       your code will render all the scenes!

const (
    scene_filenames = [
        'P02_00_sphere',
        'P02_01_sphere_ambient',
        'P02_02_sphere_room',
        'P02_03_quad',
        'P02_04_quad_room',
        'P02_05_ball_on_plane',
        'P02_06_balls_on_plane',
        'P02_07_reflections',
        'P02_08_antialiased',
        'P02_09_planar_circle',
        'P02_10_ball_on_planar_circle',
        'P02_11_cubes_on_plane',
        'P02_12_refractions',
        'P02_creative_artifact',
        'P02_creative_1',
    ]
    depth = 10
)


////////////////////////////////////////////////////////////////////////////////////////
// module aliasing to make code a little easier to read
// ex: replacing `gfx.Scene` with just `Scene`

type Point     = gfx.Point
type Vector    = gfx.Vector
type Direction = gfx.Direction
type Normal    = gfx.Normal
type Ray       = gfx.Ray
type Color     = gfx.Color
type Image     = gfx.Image

type Intersection = gfx.Intersection
type Surface      = gfx.Surface
type Scene        = gfx.Scene


////////////////////////////////////////////////////////////////////////////////////////
// functions to implement


fn intersect_ray_surface(surface Surface, ray Ray) Intersection {
    /*
        if surface's shape is a sphere
            if ray does not intersect sphere, return no intersection
            compute ray's t value(s) at intersection(s)
            if ray's t is not a valid (between ray's min and max), return no intersection
            return intersection information
            NOTE: a ray can intersect a sphere in 0, 1, or 2 different points!

        if surface's shape is a quad
            if ray does not intersect plane, return no intersection
            compute ray's t value at intersection with plane
            if ray's t is not a valid (between min and max), return no intersection
            if intersection is outside the quad, return no intersection
            return intersection information
    */
    if surface.shape == .sphere {
        e_c := surface.frame.o.vector_to(ray.e) // is equal to e - c
        a := ray.d.dot(ray.d) // a in the quadratic formula 
        b := 2 * ray.d.dot(e_c) // b = 2d^ * (e - c)
        c := e_c.dot(e_c) - math.pow(surface.radius, 2) // 1.0 // || e - c ||^2 - r^2 
        d := ( b * b) - (4 * a * c)
        if d < 0 {
            return gfx.no_intersection
        }

        t1 := (-b - math.sqrt(d)) / (2 * a)
        t2 := (-b + math.sqrt(d)) / (2 * a)
        
        if t1 >= ray.t_min && t1 <= ray.t_max {
             // calculation for the point to get normal 
            p := ray.e.add(ray.d.scale(t1))
            normal := surface.frame.o.direction_to(p) 

            return gfx.Intersection{
                frame : gfx.frame_oz(p, normal)
                surface : surface
                distance : t1
            }
        } else if t2 >= ray.t_min && t2 <= ray.t_max {
            p := ray.e.add(ray.d.scale(t2))
            normal := surface.frame.o.direction_to(p) 

            return gfx.Intersection{
                frame : gfx.frame_oz(p, normal)
                surface : surface
                distance : t2
            }
        }

    }if surface.shape == .quad {
        n := surface.frame.z
        if ray.d.dot(n) == 0 {
            return gfx.no_intersection
        }
        t := ((ray.e.vector_to(surface.frame.o)).dot(n)) / (ray.d.dot(n))
        p := ray.e.add(ray.d.scale(t))
        if t >= ray.t_min && t <= ray.t_max {
            l := p.vector_to(surface.frame.o).linf_norm()
            if l <= surface.radius {
                return gfx.Intersection {
                    frame : gfx.frame_oz(p, n)
                    surface : surface
                    distance : t
                }
            } else {
                return gfx.no_intersection
            }
        } else {
            return gfx.no_intersection
        }
    }

    if surface.shape == .planar_circle {
        n := surface.frame.z
        if ray.d.dot(n) == 0 {
            return gfx.no_intersection
        }
        t := ((ray.e.vector_to(surface.frame.o)).dot(n)) / (ray.d.dot(n))
        p := ray.e.add(ray.d.scale(t))
        if t >= ray.t_min && t <= ray.t_max {
            l := p.distance_to(surface.frame.o)
            if l <= surface.radius {
                return gfx.Intersection {
                    frame : gfx.frame_oz(p, n)
                    surface : surface
                    distance : t
                }
            } else {
                return gfx.no_intersection
            }
        } else {
            return gfx.no_intersection
        }
    }



    if surface.shape == .cube {
        mut closest := gfx.no_intersection
        size := surface.radius / 2
        centers := [
            surface.frame.o.add(surface.frame.x.scale(size)), // right face
            surface.frame.o.sub(surface.frame.x.scale(size)), // left face
            surface.frame.o.add(surface.frame.y.scale(size)), // top face
            surface.frame.o.sub(surface.frame.y.scale(size)), // bottom face
            surface.frame.o.add(surface.frame.z.scale(size)), // front face
            surface.frame.o.sub(surface.frame.z.scale(size)), // back face
        ]
        normals := [
            surface.frame.x,  // right face
            surface.frame.x.negate(),  // left face
            surface.frame.y,  // top face
            surface.frame.y.negate(),  // bottom face
            surface.frame.z,  // front face
            surface.frame.z.negate(),  // back face
        ]
        for i in 0 .. 6 {
            normal := normals[i]
            center := centers[i]
            if ray.d.dot(normal) == 0 {
                continue
            }
            t := ((ray.e.vector_to(center)).dot(normal)) / (ray.d.dot(normal))
            p := ray.e.add(ray.d.scale(t))
            if t >= ray.t_min && t <= ray.t_max {
                l := p.vector_to(center).linf_norm()
                if l <= size {
                    if closest.miss() || t < closest.distance {
                        closest = gfx.Intersection {
                            frame : gfx.frame_oz(p, normal)
                            surface : surface
                            distance : t
                        }
                    }
                } 
            }
        }
        return closest
    }
    return gfx.no_intersection
}

// Determines if given ray intersects any surface in the scene.
// If ray does not intersect anything, null is returned.
// Otherwise, details of first intersection are returned as an `Intersection` struct.
fn intersect_ray_scene(scene Scene, ray Ray) Intersection {
    mut closest := gfx.no_intersection  // type is Intersection

    /* 
        for each surface in surfaces
            continue if ray did not hit surface ( ex: inter.miss() )
            continue if new intersection is not closer than previous closest intersection
            set closest intersection to new intersection
    */
    for surface in scene.surfaces {
        inter := intersect_ray_surface(surface, ray)
        if inter.miss() {
            continue
        }
        if closest.miss() || inter.distance < closest.distance {
            closest = inter
        }
    } 
   
    return closest  // return closest intersection
}

// const (
//     depth = 8
// )

// Computes irradiance (as Color) from scene along ray
fn irradiance(scene Scene, ray Ray, depth int) Color {
    if depth > 10 {
        return gfx.black
    }
    mut accum := gfx.black

    /*
        get scene intersection
        if not hit, return scene's background intensity
        accumulate color starting with ambient
        foreach light
            compute light response    (L)
            compute light direction   (l)
            compute light visibility  (V)
            compute material response (BRFD*cos) and accumulate
        if material has reflections (lightness of kr > 0)
            create reflection ray
            accumulate reflected light (recursive call) scaled by material reflection
        return accumulated color
    */
    inter := intersect_ray_scene(scene, ray) // get scene intersecrtion
    if inter.miss() { // if not hit, returns scene's background intesity
        return scene.background_color
    }
    ambient := inter.surface.material.kd.mult(scene.ambient_color) //inter.surface.mult(scene.ambient_color // accumulate color startting with ambient
    accum = accum.add(ambient)
    
    // foreach light
    for light in scene.lights {
        s_p := light.frame.o.distance_squared_to(inter.frame.o)
        light_response := light.kl.scale(1 / s_p)
        light_direction := inter.frame.o.direction_to(light.frame.o) // compute light direction
        shadow_ray := gfx.Ray {
            e : inter.frame.o
            d : light_direction
            t_min : ray.t_min
            t_max : inter.frame.o.distance_to(light.frame.o)
        }
        // light visibiity 
        shadow_inter := intersect_ray_scene(scene, shadow_ray)
        if shadow_inter.hit() && shadow_inter.distance < light_direction.magnitude() {
            continue 
        } 

        mut occluded := 0
        if shadow_inter == gfx.no_intersection {
            occluded = 1
        }

        //compute material response and accumulate 
        //l := inter.frame.o.vector_to(light.frame.o) //inter.frame.o // might chnage to direction
        //v := inter.frame.o.vector_to(ray.e)
        //l_v := (l.add(v).as_direction()) // l_v := l.as_direction().add(v.as_direction())        
        //h := l.add(v).normalize()
        l := inter.frame.o.direction_to(light.frame.o) // l as direction
        v := inter.frame.o.direction_to(ray.e) // v as direction
        //h := l.add(v).normalize() // halfway vector
        h := (l.as_vector() + v.as_vector()).as_direction()
        // material response
        brdf := (inter.surface.material.ks.scale(math.pow(math.max(0.0, inter.frame.z.dot(h)), inter.surface.material.n)))

        accum = accum + light_response.mult(inter.surface.material.kd.add(brdf)).scale(math.abs(inter.frame.z.dot(l)) * occluded)

    } 

    if  inter.surface.material.kr.lightness() > 0 {
        reflect_direction := ray.d.negate().reflect(inter.frame.z)
        reflect_ray := gfx.Ray {
            e : inter.frame.o
            d : reflect_direction
            t_min : ray.t_min //f64(10e-5)
            t_max : ray.t_max // f64(math.inf(1))
        }
        reflect_color := irradiance(scene, reflect_ray, 0)
        accum = accum + reflect_color.mult(inter.surface.material.kr)

    }    
    return accum
}


// Computes image of scene using basic Whitted raytracer.
fn raytrace(scene Scene) Image {
    mut image := gfx.Image.new(scene.camera.sensor.resolution)

    /*
        if no anti-aliasing
            foreach image row (scene.resolution.height)
                foreach pixel in row (scene.resolution.width)
                    compute ray-camera parameters (u,v) for pixel
                    compute camera ray through pixel
                    intersect ray with scene
                    set pixel to color raytraced with ray (`irradiance`)
        else
            foreach image row
                foreach pixel in row
                    init accumulated color
                    foreach sample in y
                        foreach sample in x
                            compute ray-camera parameters (u,v) for pixel sample
                            compute camera ray through pixel sample
                            intersect ray with scene
                            accumulate color raytraced with ray (`irradiance`)
                    set pixel to average of accum color (scale by total number of samples)
        return rendered image
    */
    if scene.camera.sensor.samples == 1 {
        for y in 0 .. scene.camera.sensor.resolution.height{
            for x in 0 .. scene.camera.sensor.resolution.width {
                u := f64 (x) / scene.camera.sensor.resolution.width
                v := 1 - f64 (y) / scene.camera.sensor.resolution.height
                q := scene.camera.frame.o.add(scene.camera.frame.x.scale((u - 0.5)*scene.camera.sensor.size.width)).add(scene.camera.frame.y.scale((v - 0.5)*scene.camera.sensor.size.height)).sub((scene.camera.frame.z).scale(scene.camera.sensor.distance))
                ray := scene.camera.frame.o.ray_through(q)
                camera := Ray {
                    e : ray.e
                    d : ray.d
                    t_min : f64(10e-5)
                    t_max : f64(math.inf(1))
                    
                }
                image.set_xy(x, y, irradiance(scene, camera, 0))
            }
        } 
    } else {
        scale := 1.0 / f64(scene.camera.sensor.samples * scene.camera.sensor.samples)
        for y in 0 .. scene.camera.sensor.resolution.height { 
            for x in 0 .. scene.camera.sensor.resolution.width {
                mut color_accum := gfx.black // init accumulated color 
                // for each sample in y anf for each sample in x 
                for sy in 0 .. scene.camera.sensor.samples {
                    for sx in 0 .. scene.camera.sensor.samples {
                        u := (f64 (x) + (f64 (sx) + 0.5) / f64(scene.camera.sensor.samples)) / f64(scene.camera.sensor.resolution.width) //(scene.camera.sensor.samples) 
                        v := 1.0 - (f64 (y) + (f64 (sy) + 0.5) / f64(scene.camera.sensor.samples)) / f64(scene.camera.sensor.resolution.height) //(scene.camera.sensor.samples)
                        q := scene.camera.frame.o.add(scene.camera.frame.x.scale((u - 0.5)*scene.camera.sensor.size.width)).add(scene.camera.frame.y.scale((v - 0.5)*scene.camera.sensor.size.height)).sub((scene.camera.frame.z).scale(scene.camera.sensor.distance))
                        ray := scene.camera.frame.o.ray_through(q)
                        camera := Ray {
                            e : ray.e
                            d : ray.d
                            t_min : f64(10e-5)
                            t_max : f64(math.inf(1))
                    
                        }
                        color_accum = color_accum + irradiance(scene, camera, 0)
                    }
                } 
                image.set_xy(x, y, color_accum.scale(scale))
            }
        }
    }

    return image
}


fn main() {
    // Make sure images folder exists, because this is where all generated images will be saved
    if !os.exists('output') {
        os.mkdir('output') or { panic(err) }
    }

    for filename in scene_filenames {
        println('Rendering ${filename}...')
        scene := gfx.scene_from_file('scenes/${filename}.json')!
        image := raytrace(scene)
        image.save_png('output/${filename}.png')
    }

    println('Done!')
}
