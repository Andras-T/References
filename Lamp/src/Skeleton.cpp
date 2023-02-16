#include "framework.h"

template<class T> struct Dnum {
    float f;
    T d;
    Dnum(float f0 = 0, T d0 = T(0)) { f = f0, d = d0; }
    Dnum operator+(Dnum r) { return Dnum(f + r.f, d + r.d); }
    Dnum operator-(Dnum r) { return Dnum(f - r.f, d - r.d); }
    Dnum operator*(Dnum r) {
        return Dnum(f * r.f, f * r.d + d * r.f);
    }
    Dnum operator/(Dnum r) {
        return Dnum(f / r.f, (r.f * d - r.d * f) / r.f / r.f);
    }
};

template<class T> Dnum<T> Exp(Dnum<T> g) { return Dnum<T>(expf(g.f), expf(g.f)*g.d); }
template<class T> Dnum<T> Sin(Dnum<T> g) { return  Dnum<T>(sinf(g.f), cosf(g.f)*g.d); }
template<class T> Dnum<T> Cos(Dnum<T>  g) { return  Dnum<T>(cosf(g.f), -sinf(g.f)*g.d); }
template<class T> Dnum<T> Tan(Dnum<T>  g) { return Sin(g) / Cos(g); }
template<class T> Dnum<T> Sinh(Dnum<T> g) { return  Dnum<T>(sinh(g.f), cosh(g.f)*g.d); }
template<class T> Dnum<T> Cosh(Dnum<T> g) { return  Dnum<T>(cosh(g.f), sinh(g.f)*g.d); }
template<class T> Dnum<T> Tanh(Dnum<T> g) { return Sinh(g) / Cosh(g); }
template<class T> Dnum<T> Log(Dnum<T> g) { return  Dnum<T>(logf(g.f), g.d / g.f); }
template<class T> Dnum<T> Pow(Dnum<T> g, float n) {
    return  Dnum<T>(powf(g.f, n), n * powf(g.f, n - 1) * g.d);
}

typedef Dnum<vec2> Dnum2;

const int tessellationLevel = 60;

struct Camera {
    vec3 wEye, wLookat, wVup;
    float fov, asp, fp, bp;
public:
    Camera() {
        asp = (float)windowWidth / windowHeight;
        fov = 75.0f * (float)M_PI / 180.0f;
        fp = 5; bp = 100;
    }
    mat4 V() {
        vec3 w = normalize(wEye - wLookat);
        vec3 u = normalize(cross(wVup, w));
        vec3 v = cross(w, u);
        return TranslateMatrix(wEye * (-1)) * mat4(u.x, v.x, w.x, 0,
                                                   u.y, v.y, w.y, 0,
                                                   u.z, v.z, w.z, 0,
                                                   0,   0,   0,   1);
    }

    mat4 P() {
        return mat4(1 / (tan(fov / 2)*asp), 0,                0,                      0,
                    0,                      1 / tan(fov / 2), 0,                      0,
                    0,                      0,                -(fp + bp) / (bp - fp), -1,
                    0,                      0,                -2 * fp*bp / (bp - fp),  0);
    }
};

struct Material {
    vec3 kd, ks, ka;
    float shininess;
};

struct Light {
    vec3 La, Le;
    vec4 wLightPos;
};

struct RenderState {
    mat4	           MVP, M, Minv, V, P;
    Material *         material;
    std::vector<Light> lights;
    vec3	           wEye;
};

class Shader : public GPUProgram {
public:
    virtual void Bind(RenderState state) = 0;

    void setUniformMaterial(const Material& material, const std::string& name) {
        setUniform(material.kd, name + ".kd");
        setUniform(material.ks, name + ".ks");
        setUniform(material.ka, name + ".ka");
        setUniform(material.shininess, name + ".shininess");
    }

    void setUniformLight(const Light& light, const std::string& name) {
        setUniform(light.La, name + ".La");
        setUniform(light.Le, name + ".Le");
        setUniform(light.wLightPos, name + ".wLightPos");
    }
};

class PhongShader : public Shader {
    const char * vertexSource = R"(
		#version 330
		precision highp float;

		struct Light {
			vec3 La, Le;
			vec4 wLightPos;
		};

		uniform mat4  MVP, M, Minv;
		uniform Light[8] lights;
		uniform int   nLights;
		uniform vec3  wEye;

		layout(location = 0) in vec3  vtxPos;
		layout(location = 1) in vec3  vtxNorm;
		layout(location = 2) in vec2  vtxUV;

		out vec3 wNormal;
		out vec3 wView;
		out vec3 wLight[8];
		out vec2 texcoord;
        out vec3 vtxNoraml;

		void main() {
			gl_Position = vec4(vtxPos, 1) * MVP;
			vec4 wPos = vec4(vtxPos, 1) * M;
			for(int i = 0; i < nLights; i++) {
				wLight[i] = lights[i].wLightPos.xyz * wPos.w - wPos.xyz * lights[i].wLightPos.w;
			}
		    wView  = wEye * wPos.w - wPos.xyz;
		    wNormal = (Minv * vec4(vtxNorm, 0)).xyz;
		    texcoord = vtxUV;
            vtxNoraml = vtxNorm;
		}
	)";

    const char * fragmentSource = R"(
		#version 330
		precision highp float;

		struct Light {
			vec3 La, Le;
			vec4 wLightPos;
		};

		struct Material {
			vec3 kd, ks, ka;
			float shininess;
		};

		uniform Material material;
		uniform Light[8] lights;
		uniform int   nLights;

		in  vec3 wNormal;
		in  vec3 wView;
		in  vec3 wLight[8];
		in  vec2 texcoord;
        in  vec3 vtxNoraml;

        out vec4 fragmentColor;

		void main() {
			vec3 N = normalize(wNormal);
			vec3 V = normalize(wView);
			if (dot(N, V) < 0) N = -N;

			vec3 ka = material.ka;
			vec3 kd = material.kd;

			vec3 radiance = vec3(0, 0, 0);
			for(int i = 0; i < nLights; i++) {
				vec3 L = normalize(wLight[i]);
				vec3 H = normalize(L + V);
				float cost = max(dot(N,L), 0), cosd = max(dot(N,H), 0);
				radiance += ka * lights[i].La +
                           (kd * cost + material.ks * pow(cosd, material.shininess)) * lights[i].Le;
			}
			fragmentColor = vec4(radiance, 1);
		}
	)";
public:
    PhongShader() { create(vertexSource, fragmentSource, "fragmentColor"); }

    void Bind(RenderState state) {
        Use();
        setUniform(state.MVP, "MVP");
        setUniform(state.M, "M");
        setUniform(state.Minv, "Minv");
        setUniform(state.wEye, "wEye");
        setUniformMaterial(*state.material, "material");

        setUniform((int)state.lights.size(), "nLights");
        for (unsigned int i = 0; i < state.lights.size(); i++) {
            setUniformLight(state.lights[i], std::string("lights[") + std::to_string(i) + std::string("]"));
        }
    }
};

class Geometry {
protected:
    unsigned int vao, vbo;
public:
    Geometry() {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
    }
    virtual void Draw() = 0;
    ~Geometry() {
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }
};

class ParamSurface : public Geometry {
    struct VertexData {
        vec3 position, normal;
        vec2 texcoord;
    };

    unsigned int nVtxPerStrip, nStrips;
public:
    ParamSurface() { nVtxPerStrip = nStrips = 0; }

    virtual void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) = 0;

    VertexData GenVertexData(float u, float v) {
        VertexData vtxData;
        vtxData.texcoord = vec2(u, v);
        Dnum2 X, Y, Z;
        Dnum2 U(u, vec2(1, 0)), V(v, vec2(0, 1));
        eval(U, V, X, Y, Z);
        vtxData.position = vec3(X.f, Y.f, Z.f);
        vec3 drdU(X.d.x, Y.d.x, Z.d.x), drdV(X.d.y, Y.d.y, Z.d.y);
        vtxData.normal = cross(drdU, drdV);
        return vtxData;
    }

    void create(int N = tessellationLevel, int M = tessellationLevel) {
        nVtxPerStrip = (M + 1) * 2;
        nStrips = N;
        std::vector<VertexData> vtxData;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= M; j++) {
                vtxData.push_back(GenVertexData((float)j / M, (float)i / N));
                vtxData.push_back(GenVertexData((float)j / M, (float)(i + 1) / N));
            }
        }
        glBufferData(GL_ARRAY_BUFFER, nVtxPerStrip * nStrips * sizeof(VertexData), &vtxData[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)offsetof(VertexData, position));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)offsetof(VertexData, normal));
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)offsetof(VertexData, texcoord));
    }

    void Draw() {
        glBindVertexArray(vao);
        for (unsigned int i = 0; i < nStrips; i++) glDrawArrays(GL_TRIANGLE_STRIP, i *  nVtxPerStrip, nVtxPerStrip);
    }
};

class Sphere : public ParamSurface {
public:
    Sphere() { create(); }
    void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) {
        U = U * 2.0f * (float)M_PI, V = V * (float)M_PI;
        X = Cos(U) * Sin(V); Y = Sin(U) * Sin(V); Z = Cos(V);
    }
};

class Paraboloid : public ParamSurface {
public:
    Paraboloid() { create(); }
    void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) {
        U = U * 2.0f * M_PI;
        V = V * 2.0f * M_PI;
        X = U * Cos(V); Y = U * Sin(V); Z = X*X+Y*Y + 1;
    }
};

class Plane : public ParamSurface {
public:
    Plane() { create(); }
    void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) {
        const float size = 10;
        U = U * size;
        V = V * size;
        X = U - size/2;
        Y = V - size/2;
        Z = 0.01f;
    }
};

class Cylinder : public ParamSurface {
public:
    Cylinder() { create(); }
    void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) {
        U = U * 2.0f * M_PI, V = V * 2;
        X = Cos(U);
        Y = Sin(U);
        Z = V;
    }
};

class Circle : public ParamSurface {
public:
    Circle() { create(); }
    void eval(Dnum2& U, Dnum2& V, Dnum2& X, Dnum2& Y, Dnum2& Z) {
        U = U * 2.0f * M_PI;
        V = V * 2.0f * M_PI;
        X = U * Cos(V);
        Y = U * Sin(V);
        Z = 0.01;
    }
};

struct Shape {
    Shader *   shader;
    Material * material;
    Geometry * geometry;
    vec3 scale, translation, rotationAxis, translation2, rotationAxis2;
    float rotationAngle, rotationAngle2;
    std::vector<Shape *> children;
public:
    Shape(Shader * _shader, Material * _material, Geometry * _geometry) :
            scale(vec3(1, 1, 1)), translation(vec3(0, 0, 0)), rotationAxis(0, 0, 1), rotationAngle(0),
            translation2(vec3(0, 0, 0)), rotationAxis2(0, 0, 1), rotationAngle2(0){
        shader = _shader;
        material = _material;
        geometry = _geometry;
    }

    virtual void SetModelingTransform(mat4& M, mat4& Minv, RenderState state) {
        M =  ScaleMatrix(scale) * RotationMatrix(rotationAngle2, rotationAxis2) * RotationMatrix(rotationAngle, rotationAxis) * TranslateMatrix(translation) * state.M;
        Minv =  state.Minv * TranslateMatrix(-translation) * RotationMatrix(-rotationAngle, rotationAxis) *
        RotationMatrix(-rotationAngle2, rotationAxis2) * ScaleMatrix(vec3(1 / scale.x, 1 / scale.y, 1 / scale.z));
    }

    void Draw(RenderState state) {
        mat4 M, Minv;
        SetModelingTransform(M, Minv, state);
        state.M = M;
        state.Minv = Minv;
        state.MVP = state.M * state.V * state.P;
        state.material = material;
        shader->Bind(state);
        geometry->Draw();
        for (Shape * child : children) child->Draw(state);
    }

    virtual void Animate(float tstart, float tend) {
        rotationAngle2 += tend - tstart;
    }
};

class Scene {
    std::vector<Shape *> objects;
    Camera camera;
    std::vector<Light> lights;
    Shape * sphereObject2;
    Shape * sphereObject3;
public:
    void Build() {
        Shader * phongShader = new PhongShader();
        Material * materialLamp = new Material;
        materialLamp->kd = vec3(0.3, 0.393, 0.486);
        materialLamp->kd = vec3(0.3, 0.393, 0.486);
        materialLamp->kd = vec3(0.3, 0.393, 0.486);
        materialLamp->shininess = 60;
        Material * materialPlane = new Material;
        materialPlane->kd = vec3(0.85, 0.45, 0.36);
        materialPlane->kd = vec3(0.85, 0.45, 0.36);
        materialPlane->kd = vec3(0.85, 0.45, 0.36);
        materialPlane->shininess = 40;

        Geometry * sphere = new Sphere();
        Geometry * cylinder = new Cylinder();
        Geometry * plane = new Plane();
        Geometry * paraboloid = new Paraboloid();
        Geometry *circle = new Circle();

        Shape * paraboloidObject = new Shape(phongShader, materialLamp, paraboloid);
        paraboloidObject->translation = vec3(0.0, 0, 0.0);
        paraboloidObject->scale = 0.4 * vec3(1.7f, 1.7f, 0.45f);
        paraboloidObject->rotationAngle = -M_PI/2;
        paraboloidObject->rotationAxis = vec3 (1, 0, 0);

        sphereObject3 = new Shape(phongShader, materialLamp, sphere);
        sphereObject3->translation = vec3(0.0f, 8.5, 0);
        sphereObject3->scale = 0.8 * vec3(1.0f, 1.0f, 1.0f);
        sphereObject3->rotationAxis2 = vec3(2, 10, 5);
        sphereObject3->rotationAngle = M_PI/2;
        sphereObject3->children.push_back(paraboloidObject);

        Shape * cylinderObject3 = new Shape(phongShader, materialLamp, cylinder);
        cylinderObject3->translation = vec3(0.0, 0.0, 0.0);
        cylinderObject3->scale = (1/0.5) * vec3(0.25f, 0.25f, 2.0f);
        cylinderObject3->rotationAngle = -M_PI/2;
        cylinderObject3->rotationAxis = vec3 (1, 0, 0);

        sphereObject2 = new Shape(phongShader, materialLamp, sphere);
        sphereObject2->translation = vec3(0.0f, 10, 0);
        sphereObject2->scale = vec3(1.0f, 1.0f, 1.0f);
        sphereObject2->rotationAxis2 = vec3(8, 20, 8);
        sphereObject2->children.push_back(cylinderObject3);
        sphereObject2->children.push_back(sphereObject3);

        Shape * cylinderObject2 = new Shape(phongShader, materialLamp, cylinder);
        cylinderObject2->translation = vec3(0.0, 0.75, 0.0);
        cylinderObject2->scale = (1/0.45) * vec3(0.25f, 0.25f, 2.0f);
        cylinderObject2->rotationAngle = -M_PI/2;
        cylinderObject2->rotationAxis = vec3 (1, 0, 0);

        Shape * sphereObject1 = new Shape(phongShader, materialLamp, sphere);
        sphereObject1->translation = vec3(0.0, 0.5, 0.0);
        sphereObject1->scale = 0.4 * vec3(1.0f, 1.0f, 1.0f);
        sphereObject1->rotationAxis2 = vec3(10, 30, 10);
        sphereObject1->rotationAngle = 0;
        sphereObject1->children.push_back(cylinderObject2);
        sphereObject1->children.push_back(sphereObject2);

        Shape * circleObject1 = new Shape(phongShader, materialLamp, circle);
        circleObject1->translation = vec3(0.0, 0.40, 0.0);
        circleObject1->rotationAxis = vec3(1, 0, 0);
        circleObject1->rotationAngle= M_PI/2;
        circleObject1->scale = 0.318 * vec3(1.0f, 1.0f, 1.0f);

        Shape * cylinderObject1 = new Shape(phongShader, materialLamp, cylinder);
        cylinderObject1->translation = vec3(0, 0.20, 0.0);
        cylinderObject1->scale = vec3(2.0f, 2.0f, 0.1f);
        cylinderObject1->rotationAngle = -M_PI/2;
        cylinderObject1->rotationAxis = vec3 (1, 0, 0);

        Shape * planeObject = new Shape(phongShader, materialPlane, plane);
        planeObject->translation = vec3(0.0, 0, 0.0);
        planeObject->scale = vec3(100.0f, 100.0f, 100.0f);
        planeObject->rotationAngle = M_PI/2;
        planeObject->rotationAxis = vec3 (1, 0, 0);

        camera.wEye = vec3(-15, 15, 0);
        camera.wLookat = vec3(0, 0, 0.0);
        camera.wVup = vec3(0, 1, 0);

        lights.resize(2);
        lights[0].wLightPos = vec4(5, 5, 4, 0);
        lights[0].La = 1 * vec3(1.0f, 1.0f, 1.0f);
        lights[0].Le = 1 * vec3(1, 1, 1);
        lights[1].wLightPos = vec4(-1.0, 8.0, 0.0, 1);
        lights[1].La = 1 * vec3(1.0f, 1.0f, 1.0f);
        lights[1].Le = 1.5 * vec3(1, 1, 1);
        objects.push_back(planeObject);
        objects.push_back(cylinderObject1);
        objects.push_back(circleObject1);
        objects.push_back(sphereObject1);
    }
    RenderState stateL;
    void Render() {
        RenderState state;
        state.M =  mat4(vec4(1, 0, 0, 0),
                        vec4(0, 1, 0, 0),
                        vec4(0, 0, 1, 0),
                        vec4(0, 0, 0, 1));
        state.Minv = state.M;
        state.wEye = camera.wEye;
        state.V = camera.V();
        state.P = camera.P();
        state.lights = lights;
        for (Shape * obj : objects) obj->Draw(state);
        stateL.M = state.M;
    }

    void Animate(float tstart, float tend) {
        for (Shape * obj : objects){
            obj->Animate(tstart, tend);
        }
        sphereObject2->Animate(tstart, tend);
        sphereObject3->Animate(tstart, tend);
    }
};

Scene scene;

void onInitialization() {
    glViewport(0, 0, windowWidth, windowHeight);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    scene.Build();
}

void onDisplay() {
    glClearColor(0.02f, 0.02f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    scene.Render();
    glutSwapBuffers();
}

void onKeyboard(unsigned char key, int pX, int pY) { }

void onKeyboardUp(unsigned char key, int pX, int pY) { }

void onMouse(int button, int state, int pX, int pY) { }

void onMouseMotion(int pX, int pY) {
}

void onIdle() {
    static float tend = 0;
    const float dt = 0.1f;
    float tstart = tend;
    tend = glutGet(GLUT_ELAPSED_TIME) / 600.0f;
    for (float t = tstart; t < tend; t += dt) {
        float Dt = fmin(dt, tend - t);
        scene.Animate(t, t + Dt);
    }
    glutPostRedisplay();
}