#version 330 core

out vec4 FragColor;

uniform float act;
uniform float time;
uniform vec2 resolution;

uniform vec3 cameraPosition;
uniform vec3 cameraLookAt;

uniform float rayMaxSteps;
uniform float rayThreshold;

uniform vec3 lightPosition;

uniform vec3 diskPosition;
uniform vec3 diskRotation;
uniform vec3 diskIterationTranslate;
uniform vec3 diskIterationRotation;
uniform float diskScale;
uniform float diskOffset;
uniform vec3 diskFold;
uniform float diskIterations;
uniform float diskTime;

uniform vec3 mandlePosition;
uniform vec3 mandleRotation;
uniform float mandleSize;
uniform float mandleScale;
uniform float mandleMinRad;
uniform float mandleLimit;
uniform float mandleFactor;
uniform float mandleIterations;
uniform float mandleFoldingLimit;
uniform float mandleRadClamp1;
uniform float mandleRadClamp2;

uniform sampler2D texture01;
uniform sampler2D texture01nm;
uniform sampler2D texture02;
uniform sampler2D texture02nm;
uniform sampler2D texture03nm;
uniform sampler2D labelTexture;
uniform sampler2D ttDiskTexture;
uniform sampler2D bootTexture;

in float[12] sines;
in float[12] coses;

#define PI 3.14159265359
struct textureOptions {
    int index;
    vec2 offset;
    vec2 scale;
    bool normalMap;
};

struct material {
    vec3 ambient;
    float ambientStrength;

    vec3 diffuse;
    float diffuseStrength;

    vec3 specular;
    float specularStrength;
    float shininess;

    float shadowHardness;
    bool receiveShadows;

    textureOptions textureOptions;
};

struct entity {
    float dist;
    vec3 point;
    bool needNormals;

    material material;
};

struct hit {
    vec3 point;
    vec3 normal;

    float steps;
    float dist;
    
    float shadow;
    float last;

    entity entity;
};



//Source http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
}

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h);
}

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h);
}

vec3 opTwist(vec3 p, float angle)
{
    float c = cos(angle * p.y);
    float s = sin(angle * p.y);
    mat2 m = mat2(c, -s, s, c);
    vec3 q = vec3(m * p.xz, p.y);
    return q;
}

vec3 opBend(vec3 p, float angle)
{
    float c = cos(angle * p.y);
    float s = sin(angle * p.y);
    mat2 m = mat2(c, -s, s, c);
    vec3 q = vec3(m * p.xy, p.z);
    return q;
}

float opRound(float p, float rad)
{
    return p - rad;
}

float opOnion(float p, float thickness) {
    return abs(p) - thickness;
}

//Distance functions to creat primitives to 3D world
//Source http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdPlane(vec3 p, vec3 pos, vec4 n)
{
  // n must be normalized
    vec3 p1 = vec3(p) + pos;
    return dot(p1, n.xyz) + n.w;
}

float sdSphere(vec3 p, vec3 pos, float radius)
{
    vec3 p1 = vec3(p) + pos;
    return length(p1) - radius;
}

float sdEllipsoid(vec3 p, vec3 pos, vec3 r)
{
    vec3 p1 = p + pos;
    float k0 = length(p1 / r);
    float k1 = length(p1 / (r * r));
    return k0 * (k0 - 1.0) / k1;
}

float sdBox(vec3 p, vec3 pos, vec3 b, float r)
{   
    vec3 p1 = vec3(p) + pos;
    vec3 d = abs(p1) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0)) - r;
}

float sdTorus(vec3 p, vec3 pos, vec2 t)
{   
    vec3 p1 = vec3(p) + pos;
    vec2 q = vec2(length(p1.xz)-t.x,p1.y);
    return length(q)-t.y;
}

float sdCylinder(vec3 p, vec3 pos, vec3 c )
{
    vec3 p1 = p + pos;
    return length(p1.xz - c.xy) - c.z;
}

float sdRoundCone(in vec3 p, vec3 pos,in float r1, float r2, float h)
{    
    vec3 p1 = vec3(p) + pos;
    vec2 q = vec2( length(p1.xz), p1.y );
    
    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q,vec2(-b,a));
    
    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q-vec2(0.0,h)) - r2;
        
    return dot(q, vec2(a,b) ) - r1;
}

float sdCapsule(vec3 p, vec3 pos, vec3 a, vec3 b, float r)
{   
    vec3 p1 = vec3(p) + pos;
    vec3 pa = p1 - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}

float sdHexPrism(vec3 p, vec3 pos, vec2 h)
{
    vec3 p1 = p + pos;
    vec3 q = abs(p1);
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}


entity mMandleBox(vec3 path, material material, float size, float scale, float minrad, float limit, float factor, int iterations, float foldingLimit, float radClamp1, float radClamp2)
{
    vec4 scalev = vec4(size) / minrad;
    float absScalem1 = abs(scale - 1.0);
    float absScaleRaisedTo1mIters = pow(abs(scale), float(1 - iterations));
    vec4 p = vec4(path, 1.0), p0 = p;
 
    for (int i = 0; i < iterations; i++)
    {
        p.xyz = clamp(p.xyz, -limit, limit) * factor - p.xyz;
        float r2 = dot(p.xyz, p.xyz);
        p *= clamp(max(minrad / r2, minrad), radClamp1, radClamp2);
        p = p * scalev + p0;
        if (r2 > foldingLimit) {
            break;
        } 
   }
   entity e;
   e.dist =  ((length(p.xyz) - absScalem1) / p.w - absScaleRaisedTo1mIters);
   e.material = material;
   e.point = p.xyz;
   return e;
}

float sdMandlebulb(vec3 p, vec3 pos, float pwr, float dis, float bail, int it) {
    vec3 z = p + pos;
 
    float dr = 1.0;
    float r = 0.0;
    float power = pwr + dis;
    for (int i = 0; i < it; i++) {
        r = length(z);
        if (r > bail) break;
        
        // convert to polar coordinates
        float theta = acos(z.z/r);
        float phi = atan(z.y,z.x);
        dr =  pow(r, power - 1.0) * power * dr + 1.0;
        
        // scale and rotate the point
        float zr = pow(r, power);
        theta = theta*power;
        phi = phi*power;
        
        // convert back to cartesian coordinates
        z = zr * vec3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
        
        z += (p + pos);
    }
    return (0.5 * log(r) * r / dr);
}

float displacement(vec3 p, vec3 m)
{
    return sin(m.x*p.x)*sin(m.y*p.y)*sin(m.z*p.z);
}

float impulse(float x, float k)
{
    float h = k * x;
    return h * exp(1.0 - h);
}

float sinc(float x, float k)
{
    float a = PI * k * x - 1.0;
    return sin(a) / a;
}

float hash(vec3 p)  // replace this by something better
{
    p  = 50.0*fract( p*0.3183099 + vec3(0.71,0.113,0.419));
    return -1.0+2.0*fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

float fbm(vec3 p, int octaves) {
    vec3 p1 = p;
    float h = 0.0, a = 1.0;    
    for (int i = 0; i < octaves; ++i) {
        h += 1.0-abs(a * hash(p1)); // ridged perlin noise
        a *= 0.45; p1 *= 2.02;
    }        
    return h;
}

float sinusoidalPlasma(in vec3 p, float t, vec3 a, float c){

    return sin(p.x + t * a.x) * cos(p.y + t * a.y) * sin(p.z + t * a.z) * c;
}

vec3 rotX(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        p.x,
        c*p.y-s*p.z,
        s*p.y+c*p.z
    );
}

vec3 rotY(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x+s*p.z,
        p.y,
        -s*p.x+c*p.z
    );
}
 

vec3 rotZ(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x-s*p.y,
        s*p.x+c*p.y,
        p.z
    );
}
vec3 rot(vec3 p, vec3 a) {
    return rotX(rotY(rotZ(p, a.z), a.y), a.x);
}
vec3 translate(vec3 p, vec3 p1) {
    return p + (p1 * -1.0);
}

vec3 scale(vec3 p, float s) {
    vec3 p1 = p;
    p1 /= s;
    return p1;
} 

vec3 repeat(vec3 p, vec3 c) {
    vec3 path1 = mod(p, c) - 0.5 * c;
    return path1;
}

entity opUnion(entity m1, entity m2) {
    return m1.dist < m2.dist ? m1 : m2;
}

entity opSubtraction(entity m1, entity m2) {
    if(-m1.dist > m2.dist) {
        m1.dist *= -1.0;
        return m1;
    }
    else {
        return m2;
    }
    
}

entity opIntersection(entity m1, entity m2) {
    return m1.dist > m2.dist ? m1 : m2;
}

vec3 planeFold(vec3 z, vec3 n, float d) {
    vec3 z1 = z;
	z1.xyz -= 2.0 * min(0.0, dot(z1.xyz, n) - d) * n;
    return z1;
}

vec3 absFold(vec3 z, vec3 c) {
    vec3 z1 = z;
	z1.xyz = abs(z1.xyz - c) + c;
    return z1;
}

vec3 sierpinskiFold(vec3 z) {
    vec3 z1 = z;
	z1.xy -= min(z1.x + z1.y, 0.0);
	z1.xz -= min(z1.x + z1.z, 0.0);
	z1.yz -= min(z1.y + z1.z, 0.0);
    return z1;
}

vec3 mengerFold(vec3 z) {
    vec3 z1 = z;
	float a = min(z1.x - z1.y, 0.0);
	z1.x -= a;
	z1.y += a;
	a = min(z1.x - z1.z, 0.0);
	z1.x -= a;
	z1.z += a;
	a = min(z1.y - z1.z, 0.0);
	z1.y -= a;
	z1.z += a;
    return z1;
}

vec3 sphereFold(vec3 z, float minR, float maxR) {
    vec3 z1 = z;
	float r2 = dot(z1.xyz, z1.xyz);
	z1 *= max(maxR / max(minR, r2), 1.0);
    return z1;
}

vec3 boxFold(vec3 z, vec3 r) {
    vec3 z1 = z;
	z1.xyz = clamp(z1.xyz, -r, r) * 2.0 - z1.xyz;
    return z1;
}

entity mCross(vec3 path, vec3 l, vec3 t, float r, float s, material material) {
    entity m;
    vec3 p1 = path;
    float d1 = sdBox(p1, vec3(0.0), vec3(t.x, t.y, l.z), r);
    float d2 = sdBox(p1, vec3(0.0), vec3(l.x, t.y, t.z), r);
    float d3 = sdBox(p1, vec3(0.0), vec3(t.x, l.y, t.z), r);
    m.dist = opSmoothUnion(d1, opSmoothUnion(d2, d3, s), s);
    m.point = p1;
    m.material = material;
    return m;
}

entity mSphere(vec3 path, float radius, material material) {
    entity m;
    vec3 p1 = path;
    m.dist = sdSphere(path, vec3(0.0), radius);
    m.point = p1;
    m.material = material;
    return m;
}

entity mBox(vec3 path, vec3 size, float r, material material) {
    entity m;
    vec3 p1 = path;
    m.dist = sdBox(path, vec3(0.0), size, r);
    m.point = p1;
    m.material = material;
    return m;
}

entity mTorus(vec3 path, vec2 dim, material material) {
    entity m;
    vec3 p1 = path;
    m.dist = sdTorus(path, vec3(0.0), dim);
    m.point = p1;
    m.material = material;
    return m;
}

entity mConsole(vec3 path, vec3 pos) {

    vec3 consoleSize = vec3(12.0, 7.3, 4.0);
    vec3 baySize = vec3(10.0, 5.0, 0.3);
    vec3 screenSize = vec3(10.5, 5.419, 0.3);
    vec3 buttonSize = vec3(0.6, 0.6, 0.6);

    material consoleMaterial = material(
        vec3(0.25, 0.25, 0.25),
        1.0,
        vec3(0.4, 0.4, 0.4),
        1.0,
        vec3(0.774597, 0.774597, 0.774597),
        1.0,
        1.4,
        1.0, 
        true,
        textureOptions(
            0,
            vec2(0.0),
            vec2(0.0),
            false
        )
    );
    material screenMaterial = material(
        vec3(1.0, 1.0, 1.0),
        1.0,
        vec3(1.0, 1.0, 1.0),
        1.0,
        vec3(0.508273, 0.508273, 0.508273),
        0.0,
        0.4,
        1.0, 
        true,
        textureOptions(
            2,
            vec2(20.0, 10.0),
            vec2(0.5, -0.1),
            false
        )
    );

    entity consoleHull = mBox(
        translate(path, pos),
        consoleSize,
        0.0,
        consoleMaterial
    );
    consoleHull.needNormals = true;

    entity bayHull = mBox(
        translate(path, translate(vec3(0.0, 5.45, -1.0), pos)),
        baySize,
        0.0,
        consoleMaterial
    );
    bayHull.needNormals = true;

    entity screenHull = mBox(
        translate(path, translate(vec3(0.0, -0.7, 3.9), pos)),
        screenSize,
        0.0,
        screenMaterial
    );
    screenHull.needNormals = true;

    entity buttonHull = mBox(
        translate(path, translate(vec3(8.0, 7.1, 1.1), pos)),
        buttonSize,
        0.0,
        consoleMaterial
    );
    buttonHull.needNormals = true;

    entity complete = opSubtraction(bayHull, opSubtraction(screenHull, opUnion(buttonHull, consoleHull)));
    return complete;

}

entity mDisk(vec3 path) {
    vec3 diskSize = vec3(8.9, 9.3, 0.3);
    vec3 holeSize = vec3(0.49, 0.381, 1.0);
    vec3 labelSize = vec3(6.943, 5.419, 0.15);
    vec3 sliderCreviceSize = vec3(6.096, 3.133, 0.3);
    vec3 sliderSize = vec3(4.741, 3.133, 0.20);
    vec3 sliderHoleSize = vec3(1.185, 2.455, 1.0);
    vec3 notchSize = vec3(1.0, 1.0, 0.4);

    material bodyMaterial = material(
        vec3(0.1, 0.1, 0.1),
        0.1,
        vec3(0.1, 0.1, 0.1),
        0.1,
        vec3(0.508273, 0.508273, 0.508273),
        0.0,
        0.4,
        1.0, 
        true,
        textureOptions(
            0,
            vec2(0.0, 0.0),
            vec2(0.0, 0.0),
            true
        )
    );

    material labelMaterial = material(
        vec3(1.0, 1.0, 1.0),
        1.0,
        vec3(1.0, 1.0, 1.0),
        0.0,
        vec3(0.508273, 0.508273, 0.508273),
        1.0,
        0.0,
        1.0, 
        true,
        textureOptions(
            1,
            vec2(15.0, 11.0),
            vec2(0.5, -0.14),
            false
        )
    );

    material sliderMaterial = material(
        vec3(0.25, 0.25, 0.25),
        1.0,
        vec3(0.4, 0.4, 0.4),
        1.0,
        vec3(0.774597, 0.774597, 0.774597),
        1.0,
        1.4,
        1.0, 
        true,
        textureOptions(
            0,
            vec2(0.0),
            vec2(0.0),
            false
        )
    );

    entity bodyHull = mBox(
        path,
        diskSize,
        0.0,
        bodyMaterial
    );
    bodyHull.needNormals = true;

    entity leftHoleHull = mBox(
        translate(path, vec3(8.128, -7.45, 0.95)),
        holeSize,
        0.0,
        bodyMaterial
    );
    leftHoleHull.needNormals = true;

    entity rightHoleHull = mBox(
        translate(path, vec3(-8.128, -7.45, 0.0)),
        holeSize,
        0.0,
        bodyMaterial
    );
    rightHoleHull.needNormals = true;

    entity labelHull = mBox(
        translate(path, vec3(0.0, -3.9, 0.3)),
        labelSize,
        0.0,
        labelMaterial
    );
    labelHull.needNormals = true;

    entity sliderCreviceHull = mBox(
        translate(path, vec3(0.8, 6.18, 0.3)),
        sliderCreviceSize,
        0.0,
        bodyMaterial
    );
    sliderCreviceHull.needNormals = true;

    entity sliderHull = opSubtraction(
        mBox(
            translate(path, vec3(-2.2, 5.8, 0.0)),
            sliderHoleSize,
            0.0,
            sliderMaterial
        ),
        mBox(
            translate(path, vec3(-0.45, 6.15, -0.1)),
            sliderSize,
            0.05,
            sliderMaterial
        )
    );
    sliderHull.needNormals = true;

    entity notchHull = mBox(
        rotZ(translate(path, vec3(-9.128, 9.45, 0.0)), 0.785398),
        notchSize,
        0.0,
        bodyMaterial
    );
    entity complete = opSubtraction(notchHull, opUnion(sliderHull, opSubtraction(sliderCreviceHull, opSubtraction(labelHull, opSubtraction(rightHoleHull, opSubtraction(leftHoleHull, bodyHull))))));
    return complete;

}

entity mFractal(vec3 path, int iter, float s, float o, material material) {
    entity m;
    vec3 p1 = path;
    float scale = s;
    float offset = o;
    for(int i = 1; i <= iter; i++) {

        //p1 = boxFold(p1, vec3(2.0, 2.0, 2.0));
        //p1 = sierpinskiFold(p1);
        p1 = sphereFold(p1, 0.01, 1.8);
        //p1 = mengerFold(p1);
        //p1 = absFold(p1, vec3(1.2, 1.2, 1.2));
        //p1 = planeFold(p1, normalize(vec3(1.0, 1.0, 1.0)), 0.5);
        //p1 = rotZ(p1, time);
        //p1 = rotY(rotZ(rotX(p1, 2.2), 0.4), 0.5);
        p1 = translate(p1, vec3(1.0, 1.0, 1.0));
        p1 *= scale - offset * (scale - 1.0);
       
    }

    m = mBox(p1, vec3(1.0, 1.0, 1.0), 0.0, material);
    //this makes further objects darker
    //m.dist *= pow(scale, -float(iter));
    return m;
}

entity mDiskFractal(vec3 path, int iter, float s, float o, vec3 fold, vec3 rotation, vec3 trans) {
    entity m;
    vec3 p1 = path;
    float scale = s;
    float offset = o;
    for(int i = 1; i <= iter; i++) {

        p1 = boxFold(p1, fold);
        
        //p1 = mengerFold(p1);
        //p1 = sphereFold(p1, -10.9, 70.0);
       
        p1 = translate(p1, trans);
        p1 = rotX(rotY(rotZ(p1, rotation.z), rotation.y), rotation.x);
       
        p1 *= scale - offset * (scale - 1.0);
       
    }
    
    m = mDisk(p1);
    //this makes further objects darker
    //m.dist *= pow(scale, -float(iter));
    m.point = p1;
    return m;
}

entity scene(vec3 path)
{   
    int a = int(act);
    if(a == 1) {
        entity disks = mDiskFractal(
            rot(translate(path, diskPosition), diskRotation),
            int(diskIterations),
            diskScale,
            diskOffset,
            diskFold,
            diskIterationRotation + sinc(diskTime, 10.0),
            diskIterationTranslate * sinc(diskTime, 4.0)
        );
        disks.needNormals = true;
        return disks;
    }
    else if(a == 2) {
        material material = material(
            vec3(0.0, 1.0, 1.0),
            1.0,
            vec3(1.0, 1.0, 1.0),
            1.0,
            vec3(0.508273, 0.508273, 0.508273),
            0.0,
            0.4,
            1.0, 
            true,
            textureOptions(
                0,
                vec2(0.0),
                vec2(0.0),
                false
            )
        );
        entity bulb;
        bulb.dist = sdMandlebulb(path, vec3(0.0), 2.0, 2.0, 2.0, 2);
        bulb.material = material;
        bulb.needNormals = true;
        
        entity mandle = mMandleBox(
            rot(translate(path, mandlePosition), mandleRotation),
            material,
            mandleSize,
            mandleScale,
            mandleMinRad, 
            mandleLimit,
            mandleFactor,
            int(mandleIterations),
            mandleFoldingLimit,
            mandleRadClamp1,
            mandleRadClamp2
        );
        mandle.needNormals = true;
        return mandle;
    }
    else if (a == 3) {
        entity console = mConsole(path, vec3(0.0));
        return console;
    }

} 

hit raymarch(vec3 rayOrigin, vec3 rayDirection) {
    vec3 eps = vec3(0.0001, 0.0, 0.0);

    hit h;
    h.steps = 0.0;
    h.last = 100.0;
    
    for(float i = 0.0; i <= rayMaxSteps; i++) {
        h.point = rayOrigin + rayDirection * h.dist;
        h.entity = scene(h.point);
        h.steps += 1.0;
        h.last = min(h.dist, h.last);
        if(h.entity.dist < rayThreshold) {
            if(h.entity.needNormals == true) {
                h.normal = normalize(vec3(
                    scene(h.point + eps.xyy).dist - scene(h.point - eps.xyy).dist,
                    scene(h.point + eps.yxy).dist - scene(h.point - eps.yxy).dist,
                    scene(h.point + eps.yyx).dist - scene(h.point - eps.yyx).dist
                ));
            }
            break;
        }
        h.dist += h.entity.dist;

    }
    
    return h;
}

vec3 ambient(vec3 color, float strength) {
    return color * strength;
} 

vec3 diffuse(vec3 normal, vec3 hit, vec3 pos, vec3 color, float strength) {
    vec3 lightDir = normalize(pos - hit);
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = diff * color * strength;
    return diffuse;
}

vec3 specular(vec3 normal, vec3 eye, vec3 hit, vec3 pos, vec3 color, float strength, float shininess) {
    vec3 lightDir = normalize(pos - hit);
    vec3 viewDir = normalize(eye - hit);
    vec3 halfwayDir = normalize(lightDir + viewDir);

    float spec = pow(max(dot(normal, halfwayDir), 0.0), shininess);
    vec3 specular = strength * spec * color;
    return specular;
} 

float shadows(vec3 ro, vec3 rd, float mint, float maxt, float k) {
    float res = 1.0;
    float ph = 1e20;
    for(float t = mint; t < maxt;)
    {
        float h = scene(ro + (rd * t)).dist;
        if(h < 0.001)
            return 0.0;
        float y = h * h / (2.0 * ph);
        float d = sqrt(h * h - y * y);
        res = min(res, k * d / max(0.0, t - y));
        ph = h;
        t += h;    
    }
    return res;
}

vec2 planarMapping(vec3 p) {
    vec3 p1 = normalize(p);
    vec2 r = vec2(0.0);
    if(abs(p1.x) == 1.0) {
        r = vec2((p1.z + 1.0) / 2.0, (p1.y + 1.0) / 2.0);
    }
    else if(abs(p1.y) == 1.0) {
        r = vec2((p1.x + 1.0) / 2.0, (p1.z + 1.0) / 2.0);
    }
    else {
        r = vec2((p1.x + 1.0) / 2.0, (p1.y + 1.0) / 2.0);
    }
    return r;
}

vec2 cylindiricalMapping(vec3 p) {
    return vec2(atan(p.y / p.x), p.z);
}

vec2 scaledMapping(vec2 t, vec2 o, vec2 s) {
    return -vec2((t.x / o.x) + s.x, (t.y / o.y) + s.y);
}

vec3 processColor(hit h, vec3 rd, vec3 eye, vec2 uv, vec3 lp)
{
    if(h.steps >= rayMaxSteps) {
        //We did not hit anything
        return vec3(0.0, 0.0, 1.0);
    }
   
    vec3 depth = vec3((1.0 - smoothstep(0.0, rayMaxSteps, float(h.steps))));
    vec3 normal =  h.normal;
    if(h.entity.material.textureOptions.index == 1) {
        depth *= texture(labelTexture, scaledMapping(h.entity.point.xy, h.entity.material.textureOptions.offset, h.entity.material.textureOptions.scale)).rgb;
    }
    else if(h.entity.material.textureOptions.index == 2) {
        depth *= texture(bootTexture, scaledMapping(h.entity.point.xy, h.entity.material.textureOptions.offset, h.entity.material.textureOptions.scale)).rgb;
    }

    vec3 ambient = ambient(h.entity.material.ambient, h.entity.material.ambientStrength);
    vec3 diffuse = diffuse(normal, h.point, lp, h.entity.material.diffuse, h.entity.material.diffuseStrength);
    vec3 specular = specular(normal, eye, h.point, lp, h.entity.material.specular, h.entity.material.specularStrength, h.entity.material.shininess);
    float shadow = 1.0;
    if(h.entity.material.receiveShadows == true) {
        shadow = shadows(h.point, normalize(lp - h.point), 0.05, 0.08, h.entity.material.shadowHardness);
    }
    vec3 lights = vec3(0.0);
    lights += ambient;
    lights += diffuse;
    lights += specular;

    vec3 result = depth;
    result *= lights;
    result *= vec3(shadow);

    float gamma = 2.2;
    vec3 correct = pow(result, vec3(1.0 / gamma));
   
    return vec3(correct);
}

vec3 drawMarching() {
    float aspectRatio = resolution.x / resolution.y;
    vec2 uv = (gl_FragCoord.xy / resolution.xy) * 2.0 - 1.0;
    
    uv.x *= aspectRatio;
    float fov = 1.0;
    vec3 camPos = vec3(cameraPosition.x, cameraPosition.y, cameraPosition.z);
    vec3 forward = normalize(cameraLookAt - camPos); 
    vec3 right = normalize(vec3(forward.z, 0.0, -forward.x));
    vec3 up = normalize(cross(forward, right)); 
    
    vec3 rd = normalize(forward + fov * uv.x * right + fov * uv.y * up);
    
    vec3 ro = vec3(camPos);
 
    hit tt = raymarch(ro, rd);
    return processColor(tt, rd, ro, uv, lightPosition); 
}

void main() {
    vec3 o = drawMarching();
    FragColor = vec4(o, 1.0);
}