#include "framework.h"

const char *const fragmentSource = R"(
    #version 330 core
    precision highp float;
    in vec2 texCoord;
    uniform float dTime;
    uniform vec3 eye;
    out vec4 outColor;
    const float pi = 3.14159265359;

    struct Material {
        vec3 ka, kd, ks;
        float  S;
    };
    struct Ray {
        vec3 start, dir;
    };
    struct Hit {
        float t;
        vec3 position, normal;
        Material material;
    };
    vec4 qInv(vec4 q) {
        return vec4(-q.xyz, q.w);
    }
    vec4 qmul(vec4 q1, vec4 q2) {
        return vec4(q1.w * q2.xyz + q2.w * q1.xyz + cross(q1.xyz, q2.xyz), q1.w * q2.w - dot(q1.xyz, q2.xyz));
    }
    vec3 qRot(vec4 q, vec3 p) {
        vec4 qInv = qInv(q);
        return qmul(qmul(q, vec4(p, 0)), qInv).xyz;
    }
    vec4 quaternion(vec3 ax, float degree) {
        return vec4(ax * sin(degree / 2), cos(degree / 2));
    }

    Hit intersect_Sphere(const Ray ray, vec3 center, float radius) {
        Hit hit;
        hit.t = -1;
        vec3 dist = ray.start - center;
        float a = dot(ray.dir, ray.dir);
        float b = 2.0 * dot(dist, ray.dir);
        float c = dot(dist, dist) - radius * radius;
        float discr = b * b - 4 * a * c;
        if (discr < 0.0) {
            return hit;
        }
        hit.t = (-b - sqrt(discr)) / 2.0 / a;
        hit.position = ray.start + ray.dir * hit.t;
        hit.normal = normalize(hit.position - center);
        hit.material.ka = vec3(0.298, 0.3921, 0.4862);
        hit.material.kd = vec3(0.298, 0.3921, 0.4862);
        return hit;
    }

    float compareHits(Hit h1, Hit h2, out vec3 normal, out vec3 color) {
        if(h1.t < 0.0 && h2.t < 0.0) {
            color = vec3(0, 0, 0);
            return -1.0;
        } else if (h2.t < 0.0) {
            color = h1.material.kd;
            normal = h1.normal;
            return h1.t;
        } else if(h1.t < 0.0) {
            color = h2.material.kd;
            normal = h2.normal;
            return h2.t;
        } else {
            if (h1.t < h2.t){
                color = h1.material.kd;
                normal = h1.normal;
                return h1.t;
            } else {
                color = h2.material.kd;
                normal = h2.normal;
                return h2.t;
            }
        }
    }

    Hit intersect_Plane(const Ray ray, vec3 point, vec3 normal){
        Hit hit;
        hit.material.ka = vec3(0.8, 0.447, 0.36);
        hit.material.kd = vec3(0.8, 0.447, 0.36);
        hit.normal = normal;
        hit.t = dot(point - ray.start, normal) / dot(ray.dir, normal);
        return hit;
    }

    Hit intersect_Cir(const Ray ray, vec3 center, vec3 normal, float radius){
        Hit hit;
        hit.material.ka = vec3(0.298, 0.3921, 0.4862);
        hit.material.kd = vec3(0.298, 0.3921, 0.4862);
        hit.t =  dot(center - ray.start, normal) / dot(ray.dir, normal);
        vec3 intersection = ray.start + ray.dir * hit.t;
        float length = length(intersection-center) * length(intersection-center);
        if (length < radius*radius){
            hit.normal = normal;
            return hit;
        }
        hit.t = -1;
        return hit;
    }

    Hit intersect_Cylinder(const Ray ray, vec3 center, float radius, float height) {
        Hit hit;
        vec2 dist = ray.start.xz - center.xz;
        float a = dot(ray.dir.xz, ray.dir.xz);
        float b = 2.0 * dot(dist, ray.dir.xz);
        float c = dot(dist, dist) - radius * radius;
        float disc = b * b - 4 * a * c;
        if (disc < 0.0) {
            Hit invalid;
            invalid.t = -1;
            return invalid;
        }
        float t1 = (-b - sqrt(disc)) / 2 / a;
        vec3 hitPosition1 = ray.start + ray.dir * t1;
        float t2 = (-b + sqrt(disc)) / 2 / a;
        vec3 hitPosition2 = ray.start + ray.dir * t2;
        if (hitPosition1.y > height + center.y ||  hitPosition1.y < center.y)
            t1 = -1.0;
        if (hitPosition2.y > height + center.y || hitPosition2.y < center.y)
            t2 = -1.0;

        if (t1 < 0.0 && t2 < 0.0) {
            hit.t = -1.0;
        }
        else if (t2 < 0.0) {
            hit.t = t1;
            hit.position = hitPosition1;
        }
        else if (t1 < 0.0) {
            hit.t = t2;
            hit.position = hitPosition2;
        }
        else {
            if (t2 < t1) {
                hit.t = t2;
                hit.position = hitPosition2;
            } else {
                hit.t = t1;
                hit.position = hitPosition1;
            }
        }
        hit.normal = hit.position - center;
        hit.normal.y = 0.0;
        hit.normal = normalize(hit.normal);
        hit.material.ka = vec3(0.298, 0.3921, 0.4862);
        hit.material.kd = vec3(0.298, 0.3921, 0.4862);
        return hit;
    }

    Hit intersect_Paraboloid(const Ray ray, vec3 center, float height){
        Hit hit;
        hit.material.ka = vec3(0.298, 0.3921, 0.4862);
        hit.material.kd = vec3(0.298, 0.3921, 0.4862);
        Ray rayT = ray;
        rayT.start = rayT.start - center;
        vec2 oc = rayT.start.xz;
        vec2 dxz = rayT.dir.xz;
        float a = dot(dxz, dxz);
        float b = 2 * dot(oc, dxz) - rayT.dir.y;
        float c = dot(oc, oc) - rayT.start.y;
        float disc = b * b - 4 * a * c;
        if (disc < 0.0){
            hit.t = -1;
            return hit;
        }
        float t1 = (-b + sqrt(disc)) / 2 / a;
        float t2 = (-b - sqrt(disc)) / 2 / a;
        vec3 hitPosition1 = ray.start + ray.dir * t1;
        vec3 hitPosition2 = ray.start + ray.dir * t2;
        if (hitPosition1.y > height + center.y)
            t1 = -1.0;
        if (hitPosition2.y > height + center.y)
            t2 = -1.0;
        if(t1 < 0 && t2 < 0){
            hit.t = -1;
            return hit;
        }else if (t1 > 0 && t2 < 0){
            hit.t = t1;
        }else if (t2 > 0 && t1 < 0){
            hit.t = t2;
        }else if (t1 < t2){
            hit.t = t1;
        }else{
            hit.t = t2;
        }
        float x = vec3(rayT.start + hit.t * rayT.dir).x;
        float z = vec3(rayT.start + hit.t * rayT.dir).z;
        hit.normal = normalize(cross(vec3(1, 2 * x, 0), vec3(0, 2 * z, 1)));
        hit.position = ray.start + ray.dir * hit.t;
        return hit;
    }
    float time = dTime / 55;
    vec3 lightPos2;

    Hit intersect(const Ray ray){
        vec4 quaternion1 = quaternion(normalize(vec3(10, 35, 10)), time);
        vec4 quaternion2 = quaternion(normalize(vec3(8, 10, -1)), time);
        vec4 quaternion3 = quaternion(normalize(vec3(2, 1.5, 2)), time);
        vec3 rO1 = qRot(quaternion1, ray.start);
        vec3 rD1 = qRot(quaternion1, ray.dir);
        vec3 rO2 = qRot(quaternion2, rO1);
        vec3 rD2 = qRot(quaternion2, rD1);
        vec3 rO3 = qRot(quaternion3, rO2);
        vec3 rD3 = qRot(quaternion3, rD2);

        vec3 cylC1 = vec3(0, 0.55, 0);
        float cylH1 = 3.5;
        float cylR1 = 0.2;
        vec3 cylC2 = vec3(0, 0, 0);
        float cylH2 = 0.3;
        float cylR2 = 2.3;
        vec3 cylC3 = qRot(quaternion2,(cylC1 + vec3(0, cylH1, 0)));
        float cylH3 = cylH1;
        float cylR3 = cylR1;

        vec3 cirC = cylC2 + vec3(0, cylH2, 0);
        vec3 cirN = vec3(0, 1, 0);
        float cirR = cylR2;

        float sphR = 0.35;
        vec3 sphPos1 = cylC1- vec3(0,sphR-0.15,0);
        vec3 sphPos2 = cylC1 + vec3(0, cylH1, 0);
        vec3 sphPos3 = cylC3 + vec3(0, cylH1, 0);

        vec3 planeNormal = vec3(0, 1, 0);
        vec3 planePoint = vec3 (0, 0, 0);

        vec3 parabCenter = qRot(quaternion3,(cylC3 + vec3(0, 0.15+cylH1, 0)));

        lightPos2 = parabCenter + vec3(0, 1, 0);
        lightPos2 = qRot(qInv(quaternion3), lightPos2);
        lightPos2 = qRot(qInv(quaternion2), lightPos2);
        lightPos2 = qRot(qInv(quaternion1), lightPos2);

        Hit parab = intersect_Paraboloid(Ray(rO3, rD3), parabCenter, 2.0);
        Hit plane = intersect_Plane(ray, planePoint, planeNormal);

        Hit sph1 = intersect_Sphere(ray, sphPos1, sphR);
        Hit sph2 = intersect_Sphere(Ray(rO1, rD1), sphPos2, sphR);
        Hit sph3 = intersect_Sphere(Ray(rO2, rD2), sphPos3, sphR);

        Hit cyl1 = intersect_Cylinder(Ray(rO1, rD1), cylC1, cylR1, cylH1);
        Hit cyl2 = intersect_Cylinder(ray, cylC2, cylR2, cylH2);
        Hit cyl3 = intersect_Cylinder(Ray(rO2, rD2), cylC3, cylR3, cylH3);

        Hit cirlcle = intersect_Cir(ray, cirC, cirN, cirR);

        sph2.normal = qRot(qInv(quaternion1), sph2.normal);
        sph3.normal = qRot(qInv(quaternion2), sph3.normal);
        sph3.normal = qRot(qInv(quaternion1), sph3.normal);

        cyl1.normal = qRot(qInv(quaternion1), cyl1.normal);
        cyl3.normal = qRot(qInv(quaternion2), cyl3.normal);
        cyl3.normal = qRot(qInv(quaternion1), cyl3.normal);

        parab.normal = qRot(qInv(quaternion3), parab.normal);
        parab.normal = qRot(qInv(quaternion2), parab.normal);
        parab.normal = qRot(qInv(quaternion1), parab.normal);

        Hit hit;
        Hit first;
        first.t = -1;

        hit.t = compareHits(first, sph1, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, sph2, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, cyl1, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, sph3,  hit.normal, hit.material.kd);
        hit.t = compareHits(hit, cyl2, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, cyl3, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, cirlcle, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, plane, hit.normal, hit.material.kd);
        hit.t = compareHits(hit, parab, hit.normal, hit.material.kd);

        hit.material.ka = hit.material.kd;
        hit.material.ks = vec3(1.0, 1.0, 1.0);

        return hit;
    }
    const float epsilon = 0.0001;
    void main() {
        vec3 lightPos = vec3(2, 15, 5);
        vec3 lookat = vec3(0,0,0);
        vec3 w = eye - lookat;
        float f = length(w);
        vec3 right = normalize(cross(vec3(0.0f, 1.0f, 0.0f), w)) * f * tan(pi / 2);
        vec3 up = normalize(cross(w, right)) * f * tan(pi / 2);
        vec3 coord = normalize(vec3(texCoord * 2 - 1, -tan(pi / 2.0)));
        vec3 p = lookat + right * coord.x + up * coord.y;

        Ray ray;
        ray.start = eye;
        ray.dir = normalize(p - eye);
        Hit hit = intersect(ray);
        if(dot(hit.normal, ray.dir) > 0.0)
            hit.normal *=-1;
        vec3 hitPosition = eye + ray.dir * hit.t;

        float dL1 = length(lightPos - hitPosition);
        float dL2 = length(lightPos2 - hitPosition);
        vec3 lC = vec3(3,3,3);
        vec3 toLight1 = (lightPos - hitPosition) / dL1;
        vec3 toLight2 = (lightPos2 - hitPosition) / dL2;

        Hit light1 = intersect(Ray(hitPosition + hit.normal * epsilon, toLight1));
        Hit light2 = intersect(Ray(hitPosition + hit.normal * epsilon, toLight2));

        float lI1 = 30.0;
        float lI2 = 45.0;
        vec3 hW1 = 2 * hit.normal * (hit.normal * toLight1) - toLight1;
        vec3 hW2 = 2 * hit.normal * (hit.normal * toLight2) - toLight2;
        float S = 250.0;

        float Ls1 = pow(max(0,dot(hW1, hit.normal)),S);
        float Ls2 = pow(max(0,dot(hW2, hit.normal)),S);
        float amb = 0.25;

        if (light1.t < dL1 && light1.t > 0.0)
            lI1 = 0.0;
        if (light2.t < dL2 && light2.t > 0.0)
            lI2 = 0.0;

        float cosTheta1 = max(dot(toLight1, hit.normal), 0.0);
        float cosTheta2 = max(dot(toLight2, hit.normal), 0.0);
        outColor = vec4((hit.material.ks * Ls1 + 1/pi * hit.material.kd) * lC * lI1 * cosTheta1 / pow(dL1, 2.0) +
        (hit.material.ks * Ls2 + 1/pi * hit.material.kd) * lC * lI2 * cosTheta2 / pow(dL2, 2.0) + amb * hit.material.ka,1);
    }
)";

const char *const vertexSource = R"(
#version 330
precision highp float;
in vec4 p;
out vec2 texCoord;

void main() {
  gl_Position = p;
  texCoord = vec2(p.x + 1, p.y + 1)/2;
}
)";

GPUProgram gpuProgram(false);
unsigned int vao;
unsigned int vbo;
int dTime = 218;
vec3 lookat(0,0,0);
vec3 eye = vec3(0, 12, 8);

void onInitialization() {
    glViewport(0, 0, windowWidth, windowHeight);
    float vertices[] = {-1.0f, -1.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f};
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    gpuProgram.create(vertexSource, fragmentSource, "outColor");
}

void onDisplay() {
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    gpuProgram.setUniform((float)dTime, "dTime");
    gpuProgram.setUniform((vec3)eye, "eye");
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glutSwapBuffers();
}

void onKeyboard(unsigned char key, int pX, int pY) {}

void onKeyboardUp(unsigned char key, int pX, int pY) {}

void onMouseMotion(int pX, int pY) {}

void onMouse(int button, int state, int pX, int pY) {}
int t = 0;
void onIdle() {
    float dt = 0.01f;
    eye = vec3((eye.x - lookat.x) * cosf(dt) + (eye.z - lookat.z) * sinf(dt) + lookat.x, eye.y,
               -(eye.x - lookat.x) * sinf(dt) + (eye.z - lookat.z) * cosf(dt) + lookat.z);
    glutPostRedisplay();
    dTime++;
    t++;
}