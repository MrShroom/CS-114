#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#define PI 3.1415926535897932384626433832795

/*
 * Thread-safe random number generator
 */

struct RNG {
    RNG() : distrb(0.0, 1.0), engines() {}

    void init(int nworkers) {
        std::random_device rd;
        engines.resize(nworkers);
        for ( int i = 0; i < nworkers; ++i )
            engines[i].seed(rd());
    }

    double operator()() {
        int id = omp_get_thread_num();
        return distrb(engines[id]);
    }

    std::uniform_real_distribution<double> distrb;
    std::vector<std::mt19937> engines;
} rng;

/*
 * Basic data types
 */

struct Vec {
    double x, y, z;

    Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }

    Vec operator+ (const Vec &b) const  { return Vec(x+b.x, y+b.y, z+b.z); }
    Vec operator- (const Vec &b) const  { return Vec(x-b.x, y-b.y, z-b.z); }
    Vec operator* (double b) const      { return Vec(x*b, y*b, z*b); }
    bool operator== (const Vec &b) const 
    {
        const double eps = 1e-4;// 2 * std::numeric_limits<double>::epsilon();
        if (std::abs(x - b.x) > eps)
            return false;
        if (std::abs(y - b.y) > eps)
            return false;
        if (std::abs(z - b.z) > eps)
            return false;
        return true;
    }

    Vec mult(const Vec &b) const        { return Vec(x*b.x, y*b.y, z*b.z); }
    Vec& normalize()                    { return *this = *this * (1.0/std::sqrt(x*x+y*y+z*z)); }
    double dot(const Vec &b) const      { return x*b.x+y*b.y+z*b.z; }
    Vec cross(const Vec&b) const        { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }
};

struct Ray {
    Vec o, d;//orgin direction
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

struct BRDF {
    virtual Vec eval(const Vec &n, const Vec &o, const Vec &i) const = 0;
    virtual void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const = 0;
};

/*
 * Utility functions
 */

inline double clamp(double x)   {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) {
    return static_cast<int>(std::pow(clamp(x), 1.0/2.2)*255+.5);
}


/*
 * Shapes
 */

struct Sphere {
    Vec p, e;           // position, emitted radiance
    double rad;         // radius
    const BRDF &brdf;   // BRDF
    
    Sphere(double rad_, Vec p_, Vec e_, const BRDF &brdf_) :
        rad(rad_), p(p_), e(e_), brdf(brdf_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b*b-op.dot(op)+rad*rad;
        if ( det<0 ) return 0; else det = sqrt(det);
        return (t = b-det)>eps ? t : ((t = b+det)>eps ? t : 0);
    }
};


/*
 * Sampling functions
 */

inline void createLocalCoord(const Vec &n, Vec &u, Vec &v, Vec &w) {
    w = n;
    u = ((std::abs(w.x)>.1 ? Vec(0, 1) : Vec(1)).cross(w)).normalize();
    v = w.cross(u);
}


void uniformRandom(const Vec &normal,  Vec &i)  {
    double z = std::sqrt(rng());
    double r = std::sqrt(1.0 - z * z);
    double phi = 2.0 * PI * rng();
    double x = r * std::cos(phi);
    double y = r * std::sin(phi);
    Vec u, v, w;
    createLocalCoord(normal, u, v, w);
    i = ((u * x) + (v *y) + (w*z)).normalize();

    
}

/*
 * BRDFs (bidirectional reflectance distribution function)
 */

// Ideal diffuse BRDF
struct DiffuseBRDF : public BRDF {
    DiffuseBRDF(Vec kd_) : kd(kd_) {}

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const {
        return kd * (1.0/PI);
    }

    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
        
        uniformRandom(n, i);
        pdf = n.dot(i)/ PI;
    }

    Vec kd;
};

// Ideal Specular BRDF
struct SpecularBRDF : public BRDF {
    
    SpecularBRDF(Vec ks_) : ks(ks_) {}

    Vec mirroredDirection(const Vec &n, const Vec &o_0) const
    {
        return (n * 2.0 * n.dot(o_0) - o_0).normalize();
    }

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const
    {
        if(o == mirroredDirection(n,i))
            return ks * (1.0 / n.dot(o));
        return Vec();
    }

    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {

        i = mirroredDirection(n,o).normalize();
        pdf = 1.0;
    }

    Vec ks;
};

/*
 * Scene configuration
 */

// Pre-defined BRDFs
const DiffuseBRDF leftWall(Vec(.75,.25,.25)),
                  rightWall(Vec(.25,.25,.75)),
                  otherWall(Vec(.75,.75,.75)),
                  blackSurf(Vec(0.0,0.0,0.0)),
                  brightSurf(Vec(0.9,0.9,0.9));
const SpecularBRDF shinySurf(Vec(0.999, 0.999, 0.999));

// Scene: list of spheres
const Sphere spheres[] = {
    Sphere(1e5,  Vec(1e5+1,40.8,81.6),   Vec(),         leftWall),   // Left
    Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec(),         rightWall),  // Right
    Sphere(1e5,  Vec(50,40.8, 1e5),      Vec(),         otherWall),  // Back
    Sphere(1e5,  Vec(50, 1e5, 81.6),     Vec(),         otherWall),  // Bottom
    Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec(),         otherWall),  // Top
    Sphere(16.5, Vec(27,16.5,47),        Vec(),         brightSurf), // Ball 1
    Sphere(16.5, Vec(73,16.5,78),        Vec(),         shinySurf),  // Ball 2
    Sphere(5.0,  Vec(50,70.0,81.6),      Vec(50,50,50), blackSurf)   // Light
};

// Camera position & direction
const Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());


/*
 * Global functions
 */

bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres)/sizeof(Sphere), d, inf = t = 1e20;
    for ( int i = int(n); i--;) if ( (d = spheres[i].intersect(r))&&d<t ) { t = d; id = i; }
    return t<inf;
}


void luminaireSample(const Vec &point, Vec &out, Vec &ny, double &pdf, Vec &light)
{
    const Sphere myLight = spheres[7];
    double z = 2 * rng() - 1.0;
    double e = rng();
    double x = std::sqrt(1.0 - z * z) *std::sin(2 * PI * e);
    double y = std::sqrt(1.0 - z * z) *std::cos(2 * PI * e);

    out = myLight.p + (Vec(x, y, z)*myLight.rad);
    pdf = 1 / (4 * PI *myLight.rad * myLight.rad);
    ny = Vec(x, y, z);
    light = myLight.e;

}

double visiblity(const Vec &x,const Vec &y)
{
    Ray r(x, (y - x).normalize());
    double t;
    int id = 0;
    intersect(r, t, id);

    Vec intersection = r.o + r.d*t;
    if (intersection == y)
           return 1.0;
    return 0.0;

};

/*
 * KEY FUNCTION: radiance estimator
 */

Vec receivedRadiance(const Ray &r, int depth, bool flag) 
{
    double t;                                   // Distance to intersection
    int id = 0;                                 // id of intersected sphere

    if ( !intersect(r, t, id) ) return Vec();   // if miss, return black
    const Sphere &obj = spheres[id];            // the hit object

    Vec x = r.o + r.d*t;                        // The intersection point
    Vec o = (Vec() - r.d).normalize();          // The outgoing direction (= -r.d)

    Vec n = (x - obj.p).normalize();            // The normal direction
    if ( n.dot(o) < 0 ) n = n*-1.0;

    Vec Le = obj.e;                             // Emitted radiance
    const BRDF &brdf = obj.brdf;                // Surface BRDF at x

    const int rrDepth = 5;                      // Depth to start RR
    const double survivalProbability = 0.9;     // Chance of surviving RR 
    Vec y_1;                                    // Point sampled from light source
    Vec ny;                                     // normal at y_1
    double pdf1;                                // pdf at of sampling y_1
    Vec Light;                                  // Emitted radiance from light source

    /*
        Direct radiance
    */
    luminaireSample(x, y_1, ny, pdf1, Light);   // sample from light
    Vec omega_1 = (y_1 - x).normalize();        // omega from x to y
    double r_squared = (x-y_1).dot(x -y_1);

    Vec directRadiance = Light.mult( brdf.eval(n, omega_1, o))
        * visiblity(x,y_1) 
        * n.dot(omega_1) 
        * (ny.dot(omega_1*-1)
        /(r_squared * pdf1));

    /*
        Indirect radiance
    */
    double p = 1.0;
    if (depth > rrDepth)
        p = survivalProbability;
    Vec indirectRadiance;
    if (rng() < p)
    {
        Vec omega_2;
        double pdf2;
        brdf.sample(n, o, omega_2, pdf2);
        Ray y(x,omega_2.normalize());
        indirectRadiance =  receivedRadiance(y, depth + 1, false)
            .mult(brdf.eval(n, o, omega_2)) 
            * (n.dot(omega_2)
            /(pdf2*p));
    }
    if(flag)//frist call we want to include radiance of object
        return  Le + directRadiance + indirectRadiance;
    return  directRadiance + indirectRadiance;
}


/*
 * Main function (do not modify)
 */

int main(int argc, char *argv[]) {
    int nworkers = omp_get_num_procs();
    omp_set_num_threads(nworkers);
    rng.init(nworkers);

    int w = 480, h = 360, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
    Vec cx = Vec(w*.5135/h), cy = (cx.cross(cam.d)).normalize()*.5135;    
    std::vector<Vec> c(w*h);

#pragma omp parallel for schedule(dynamic, 1)
    for ( int y = 0; y < h; y++ ) {
        for ( int x = 0; x < w; x++ ) {
            const int i = (h - y - 1)*w + x;

            for ( int sy = 0; sy < 2; ++sy ) {
                for ( int sx = 0; sx < 2; ++sx ) {
                    Vec r;
                    for ( int s = 0; s<samps; s++ ) {
                        double r1 = 2*rng(), dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                        double r2 = 2*rng(), dy = r2<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                        Vec d = cx*(((sx+.5 + dx)/2 + x)/w - .5) +
                            cy*(((sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        r = r + receivedRadiance(Ray(cam.o, d.normalize()), 1, true)*(1./samps);
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
                }
            }
        }
#pragma omp critical
        fprintf(stderr,"\rRendering (%d spp) %6.2f%%",samps*4,100.*y/(h-1));
    }
    fprintf(stderr, "\n");

    // Write resulting image to a PPM file
    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for ( int i = 0; i<w*h; i++ )
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    fclose(f);
    system("convert .\\image.ppm win:");
    return 0;
}
