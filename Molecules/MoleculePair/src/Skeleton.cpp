//=============================================================================================
// Mintaprogram: Zold haromszog. Ervenyes 2019. osztol.
//
// A beadott program csak ebben a fajlban lehet, a fajl 1 byte-os ASCII karaktereket tartalmazhat, BOM kihuzando.
// Tilos:
// - mast "beincludolni", illetve mas konyvtarat hasznalni
// - faljmuveleteket vegezni a printf-et kiveve
// - Mashonnan atvett programresszleteket forrasmegjeloles nelkul felhasznalni es
// - felesleges programsorokat a beadott programban hagyni!!!!!!!
// - felesleges kommenteket a beadott programba irni a forrasmegjelolest kommentjeit kiveve
// ---------------------------------------------------------------------------------------------
// A feladatot ANSI C++ nyelvu forditoprogrammal ellenorizzuk, a Visual Studio-hoz kepesti elteresekrol
// es a leggyakoribb hibakrol (pl. ideiglenes objektumot nem lehet referencia tipusnak ertekul adni)
// a hazibeado portal ad egy osszefoglalot.
// ---------------------------------------------------------------------------------------------
// A feladatmegoldasokban csak olyan OpenGL fuggvenyek hasznalhatok, amelyek az oran a feladatkiadasig elhangzottak
// A keretben nem szereplo GLUT fuggvenyek tiltottak.
//
// NYILATKOZAT
// ---------------------------------------------------------------------------------------------
// Nev    : Toth Andras
// Neptun : FKMV9M
// ---------------------------------------------------------------------------------------------
// ezennel kijelentem, hogy a feladatot magam keszitettem, es ha barmilyen segitseget igenybe vettem vagy
// mas szellemi termeket felhasznaltam, akkor a forrast es az atvett reszt kommentekben egyertelmuen jeloltem.
// A forrasmegjeloles kotelme vonatkozik az eloadas foliakat es a targy oktatoi, illetve a
// grafhazi doktor tanacsait kiveve barmilyen csatornan (szoban, irasban, Interneten, stb.) erkezo minden egyeb
// informaciora (keplet, program, algoritmus, stb.). Kijelentem, hogy a forrasmegjelolessel atvett reszeket is ertem,
// azok helyessegere matematikai bizonyitast tudok adni. Tisztaban vagyok azzal, hogy az atvett reszek nem szamitanak
// a sajat kontribucioba, igy a feladat elfogadasarol a tobbi resz mennyisege es minosege alapjan szuletik dontes.
// Tudomasul veszem, hogy a forrasmegjeloles kotelmenek megsertese eseten a hazifeladatra adhato pontokat
// negativ elojellel szamoljak el es ezzel parhuzamosan eljaras is indul velem szemben.
//=============================================================================================
#include "framework.h"

const char * const vertexShaderSource = R"(
	#version 330
	precision highp float;

    layout(location = 0) in vec2 position;
	uniform mat4 MVP;

	void main() {
		vec4 mvPosition = vec4(position.x, position.y, 0.0f, 1.0f) * MVP;
        vec4 calc = vec4 (mvPosition.x, mvPosition.y, sqrt(mvPosition.x * mvPosition.x + mvPosition.y * mvPosition.y + 1.0f), 1.0f);
		gl_Position = vec4(calc.x/(calc.z+1.0f), calc.y/(calc.z+1.0f), 0.0f, 1.0f);
	}
)";

const char * const fragmentShaderSource = R"(
	#version 330
	precision highp float;

	uniform vec3 color;
	out vec4 fragmentColor;

	void main() {
		fragmentColor = vec4(color, 1.0f);
	}
)";

GPUProgram gpuProgram;
unsigned int buffer;
unsigned int vertexArray;
unsigned int program;
long lastUpdate = 0;
float sx1 = 1.0f, sy1 = 1.0f;
float sx2 = 1.0f, sy2 = 1.0f;
float phi1 = 0.0f, phi2 = 0.0f;
vec2 wTranslate1(0.0f,0.0f), wTranslate2(0.0f,0.0f), camera2D(0,0);
const float U_CHARGE = 0.0006f;
const int CIRCLE_RESOLUTION = 32;
const float  RADIUS = 2.50f;
const float MAX_MASS = 0.03f;
const float RESISTANCE = 1.8f;
const float MIN_DST = 0.1;

int ran(int low, int high) { return low + (rand() % (high - low + 1)); }
float ranf(float low, float high) { return low + (rand() / (RAND_MAX / (high - low))); }
mat4 M1();
mat4 M2();

struct Atom{
    vec4 position;
    vec2 relativePosition;
    vec2 coulomb;
    vec3 color;
    float charge;
    float mass;

    Atom(vec2 p, float m){
        position = vec4(p.x, p.y, 0.0f, 1.0f);
        mass = m;
    }
};

struct MoleculePair{
    std::vector<vec4> buffer;
    vec2 position1, position2, ogPosition1, ogPosition2;
    std::vector<Atom> atoms1, atoms2;
    std::vector<float> edges1, edges2;
    int numberOfAtoms1, numberOfEdges1, numberOfAtoms2, numberOfEdges2;
    float mass1 = 0.0f, inertia1 = 0.0f, tau1 = 0.0f, angular_acc1 = 0.0f, angular_vel1= 0.0f;
    float mass2 = 0.0f, inertia2 = 0.0f, tau2 = 0.0f, angular_acc2 = 0.0f, angular_vel2 = 0.0f;
    vec2 acceleration1, acceleration2, velocity1, velocity2;
    mat4 first_Molecules_Matrix;
    mat4 second_Molecules_Matrix;

    MoleculePair(){
        Create();
    }

    void Create(){
        position1 = vec2(ranf(-0.6f, 0.4f), ranf(-0.6f, 0.4f));
        position2 = vec2(ranf(-0.4f, 0.6f), ranf(-0.4f, 0.6f));
        while(length(position1-position2) < 0.25){
            position1 = vec2(ranf(-0.6f, 0.4f), ranf(-0.6f, 0.4f));
            position2 = vec2(ranf(-0.4f, 0.6f), ranf(-0.4f, 0.6f));
        }
        numberOfAtoms1 = ran(2, 8);
        numberOfAtoms2 = ran(2, 8);
        mass1 = 0; mass2 = 0;
        angular_acc1 = 0;
        angular_acc2 = 0;
        angular_vel1 = 0;
        angular_vel2 = 0;
        wTranslate1 = vec2(0,0);
        wTranslate2 = vec2(0,0);
        phi1 = 0;
        phi2 = 0;
        for (int i = 0; i < numberOfAtoms1; ++i) {
            float x = ranf(-0.4f, 0.4f), y = ranf(-0.4f, 0.4f);
            while(x < MIN_DST && x > -1*MIN_DST || y < MIN_DST && y > -1*MIN_DST)
                x = ranf(-0.4f, 0.4f), y = ranf(-0.4f, 0.4f);
            vec2 pos = vec2(x+position1.x, y+position1.y);
            atoms1.push_back(Atom (vec2(pos.x, pos.y), ranf(0.005f,  MAX_MASS)));
            mass1 += atoms1[i].mass;
        }
        for (int i = 0; i < numberOfAtoms2; ++i) {
            float x = ranf(-0.4f, 0.4f), y = ranf(-0.4f, 0.4f);
            while(x < MIN_DST && x > -1*MIN_DST || y < MIN_DST && y > -1*MIN_DST)
                x = ranf(-0.4f, 0.4f), y = ranf(-0.4f, 0.4f);
            vec2 pos = vec2(x+position2.x, y+position2.y);
            atoms2.push_back(Atom (vec2(pos.x, pos.y), ranf(0.005f, MAX_MASS)));
            mass2 += atoms2[i].mass;
        }
        SetCharge();
        MoveCenter();
        SetPhysicProp();
        RandomTreeGenerator(numberOfAtoms1, edges1);
        RandomTreeGenerator(numberOfAtoms2, edges2);
        Upload();
    }

    void SetCharge(){
        float q = 0.0f;
        for (int i = 0; i < numberOfAtoms1; i++) {
            atoms1[i].charge = ranf(-1*U_CHARGE,U_CHARGE);
            while(atoms1[i].charge < 0.1*U_CHARGE && atoms1[i].charge > -0.1*U_CHARGE){
                atoms1[i].charge = ranf(-1*U_CHARGE,U_CHARGE);}
            q+=atoms1[i].charge;
        }atoms1[numberOfAtoms1-1].charge = -1 * (q-atoms1[numberOfAtoms1-1].charge);
        q=0;
        for (int i = 0; i < numberOfAtoms2; i++) {
            atoms2[i].charge = ranf(-1*U_CHARGE,U_CHARGE);
            while(atoms2[i].charge < 0.1*U_CHARGE && atoms2[i].charge > -0.1*U_CHARGE){
                atoms2[i].charge = ranf(-1*U_CHARGE,U_CHARGE);}
            q+=atoms2[i].charge;
        }atoms2[numberOfAtoms2-1].charge = -1 * (q-atoms2[numberOfAtoms2-1].charge);
        SetColor();
    }

    void SetColor(){
        for (int i = 0; i < numberOfAtoms1; ++i) {
            if (atoms1[i].charge < 0)
                atoms1[i].color = vec3 (0,0, -1/U_CHARGE*atoms1[i].charge);
            else
                atoms1[i].color = vec3(1/U_CHARGE*atoms1[i].charge,0,0);
        }
        for (int i = 0; i < numberOfAtoms2; ++i) {
            if (atoms2[i].charge < 0)
                atoms2[i].color = vec3 (0,0, -1/U_CHARGE*atoms2[i].charge);
            else
                atoms2[i].color = vec3(1/U_CHARGE*atoms2[i].charge,0,0);
        }
    }

    void MoveCenter(){
        float xShift = 0.0f, yShift = 0.0f;
        for (int i = 0; i < numberOfAtoms1; i++) {
            xShift += atoms1[i].mass * atoms1[i].position.x;
            yShift += atoms1[i].mass * atoms1[i].position.y;
        }
        position1.x = xShift * (1 / mass1);
        position1.y = yShift * (1 / mass1);
        xShift = 0; yShift = 0;
        for (int i = 0; i < numberOfAtoms2; i++) {
            xShift += atoms2[i].mass * atoms2[i].position.x;
            yShift += atoms2[i].mass * atoms2[i].position.y;
        }
        position2.x = xShift * (1 / mass2);
        position2.y =  yShift * (1 / mass2);
    }

    void SetPhysicProp(){
        for (int i = 0; i < numberOfAtoms1; ++i){
            atoms1[i].relativePosition = vec2(atoms1[i].position.x, atoms1[i].position.y) - position1;
            inertia1 += atoms1[i].mass * length(atoms1[i].relativePosition)* length(atoms1[i].relativePosition);
        }for (int i = 0; i < numberOfAtoms2; ++i){
            atoms2[i].relativePosition = vec2(atoms2[i].position.x, atoms2[i].position.y)  - position2;
            inertia2 += atoms2[i].mass * length(atoms2[i].relativePosition)* length(atoms2[i].relativePosition);
        }
    }

    void RandomTreeGenerator(int numberOfAtoms, std::vector<float> &edges){
        int length = numberOfAtoms - 2;
        int prufer[length];
        for (int i = 0; i < length; i++)
            prufer[i] = ran(0, length + 1) + 1;
        int vertex_set[numberOfAtoms];
        for (int i = 0; i < numberOfAtoms; i++)
            vertex_set[i] = 0;
        for (int i = 0; i < numberOfAtoms - 2; i++)
            vertex_set[prufer[i] - 1] += 1;
        int j = 0, k = 0;
        for (int i = 0; i < numberOfAtoms -2; i++){
            for (j = 0; j < numberOfAtoms; j++){
                if (vertex_set[j] == 0){
                    edges.push_back(j);
                    edges.push_back(prufer[i]-1);
                    vertex_set[j] = -1;
                    vertex_set[prufer[i] - 1]--;
                    break;
                }   }   }
        j = 0;
        int x, y;
        for (int i = 0; i < numberOfAtoms; i++){
            if (vertex_set[i] == 0 && j == 0){
                x = i + 1;
                j++;
            }
            else if (vertex_set[i] == 0 && j == 1)
                y = i + 1;
        }
        edges.push_back(x - 1);
        edges.push_back(y - 1);
    }

    void Calculate(){
        acceleration1 = 0; acceleration2 = 0;
        tau1 = 0; tau2 = 0;
        for (int j = 0; j <numberOfAtoms1; ++j) {
            atoms1[j].coulomb = vec2(0,0);
            for (int k = 0; k < numberOfAtoms2; ++k) {
                vec4 v = (atoms1[j].position*M1()) - (atoms2[k].position*M2());
                vec2 norm = normalize(vec2(v.x, v.y));
                vec2 tmp = (atoms1[j].charge*atoms2[k].charge) / powf(length(vec2(v.x, v.y)), 2) * norm;
                atoms1[j].coulomb.x += tmp.x;
                atoms1[j].coulomb.y += tmp.y;
            }
        }
        for (int j = 0; j <numberOfAtoms2; ++j) {
            atoms2[j].coulomb = vec2(0,0);
            for (int k = 0; k < numberOfAtoms1; ++k) {
                vec4 v = vec4(atoms2[j].position*M2()) - (atoms1[k].position*M1());
                vec2 norm = normalize(vec2(v.x, v.y));
                vec2 tmp = (atoms1[j].charge*atoms2[k].charge) / powf(length(vec2(v.x, v.y)), 2) * norm;
                atoms2[j].coulomb.x += tmp.x;
                atoms2[j].coulomb.y += tmp.y;
            }
        }
        for (int i = 0; i < numberOfAtoms1; ++i) {
            vec3 torque(cross(atoms1[i].relativePosition, atoms1[i].coulomb));
            tau1 += 1000*torque.z;
            acceleration1 = acceleration1 + atoms1[i].coulomb * (1/atoms1[i].mass);
        }
        for (int i = 0; i < numberOfAtoms2; ++i) {
            vec3 torque(cross(atoms2[i].relativePosition, atoms2[i].coulomb));
            tau2 += 1000*torque.z;
            acceleration2 = acceleration2 + atoms2[i].coulomb * (1/atoms2[i].mass);
        }
        acceleration1 = acceleration1 - velocity1 * length(velocity1) *  RESISTANCE;
        acceleration2 = acceleration2 - velocity2 * length(velocity2) * RESISTANCE;
        angular_acc1 = tau1/inertia1 - 0.01 * fabs(angular_vel1) * angular_vel1/inertia1;
        angular_acc2 = tau2/inertia2 - 0.01 * fabs(angular_vel2) * angular_vel2/inertia2;
    }

    void Upload(){
        buffer.clear();
        numberOfEdges1 = 0;
        numberOfEdges2 = 0;
        for (int i = 0; i < numberOfAtoms1-1; ++i) {
            vec4 tmp = atoms1[edges1[2*i+1]].position - atoms1[edges1[2*i]].position;
            float deltaX = (tmp.x/30), deltaY = tmp.y/30;
            for (int j = 0; j < 30; ++j) {
                vec4 tmp1 = vec4(atoms1[edges1[2*i]].position.x+deltaX*j, atoms1[edges1[2*i]].position.y+deltaY*j,0,1);
                buffer.push_back(tmp1);
                vec4 tmp2 = vec4(atoms1[edges1[2*i]].position.x+deltaX*(j+1), atoms1[edges1[2*i]].position.y+deltaY*(j+1),0,1);
                buffer.push_back(tmp2);
                numberOfEdges1++;
            }
        }
        for (int i = 0; i < numberOfAtoms2-1; ++i) {
            vec4 tmp = atoms2[edges2[2*i+1]].position - atoms2[edges2[2*i]].position;
            float deltaX = (tmp.x/20), deltaY = tmp.y/20;
            for (int j = 0; j < 20; ++j) {
                buffer.push_back(vec4(atoms2[edges2[2*i]].position.x+deltaX*j, atoms2[edges2[2*i]].position.y+deltaY*j,0,1));
                buffer.push_back(vec4(atoms2[edges2[2*i]].position.x+deltaX*(j+1), atoms2[edges2[2*i]].position.y+deltaY*(j+1),0,1));
                numberOfEdges2++;
            }
        }
        for (int j = 0; j < numberOfAtoms1; ++j) {
            float delta = M_PI * 2.0f / CIRCLE_RESOLUTION;
            for (int i = 0; i < CIRCLE_RESOLUTION; i++) {
                buffer.push_back((vec4(atoms1[j].position.x + (atoms1[j].mass*RADIUS * cosf(i * delta)),
                                       atoms1[j].position.y + (atoms1[j].mass*RADIUS * sinf(i * delta)))));
            }
        }
        for (int j = 0; j < numberOfAtoms2; ++j) {
            float delta = M_PI * 2.0f / CIRCLE_RESOLUTION;
            for (int i = 0; i < CIRCLE_RESOLUTION; i++) {
                buffer.push_back((vec4(atoms2[j].position.x + (atoms2[j].mass*RADIUS * cosf(i * delta)),
                                       atoms2[j].position.y + (atoms2[j].mass*RADIUS * sinf(i * delta)))));
            }
        }
        float delta = M_PI * 2.0f / CIRCLE_RESOLUTION;
        for (int i = 0; i < CIRCLE_RESOLUTION; i++) {
            buffer.push_back((vec4(position1.x + (0.01 * cosf(i * delta)),
                                   position1.y + (0.01 * sinf(i * delta)))));
        }
        for (int i = 0; i < CIRCLE_RESOLUTION; i++) {
            buffer.push_back((vec4(position2.x + (0.01 * cosf(i * delta)),
                                   position2.y + (0.01 * sinf(i * delta)))));
        }
        for (int i = 0; i < CIRCLE_RESOLUTION; i++) {
            buffer.push_back((vec4(0 + (0.02 * cosf(i * delta)),
                                   0 + (0.02 * sinf(i * delta)))));
        }
    }
    void Draw();
};
MoleculePair *moleculePair;
mat4 M1(){ return TranslateMatrix(vec3(-1*moleculePair->position1.x, -1*moleculePair->position1.y,0)) * RotationMatrix(phi1, vec3(0, 0, 1)) *
                  TranslateMatrix(vec3(moleculePair->position1.x, moleculePair->position1.y,0)) * TranslateMatrix(vec3(wTranslate1.x, wTranslate1.y,0)); }
mat4 M2(){ return TranslateMatrix(vec3(-1*moleculePair->position2.x, -1*moleculePair->position2.y,0)) * RotationMatrix(phi2, vec3(0, 0, 1)) *
                  TranslateMatrix(vec3(moleculePair->position2.x, moleculePair->position2.y,0)) * TranslateMatrix(vec3(wTranslate2.x, wTranslate2.y,0)); }

void onInitialization() {
    glViewport(0, 0, windowWidth, windowHeight);
    glGenVertexArrays(1,&vertexArray);
    glBindVertexArray(vertexArray);
    glGenBuffers(1, &buffer);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    moleculePair = new MoleculePair();
    glutPostRedisplay();
}

void onDisplay() { moleculePair->Draw(); }

void MoleculePair::Draw(){
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec4) * (CIRCLE_RESOLUTION * (numberOfAtoms1+numberOfAtoms2+3) + (numberOfEdges1 + numberOfEdges2)*2), &buffer[0], GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(vec4), (void*)0);
    gpuProgram.create(vertexShaderSource, fragmentShaderSource, "fragmentColor");
    glClearColor(0.42,0.42,0.42,0);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindVertexArray(vertexArray);
    gpuProgram.Use();
    glLineWidth( 2 );
    mat4 m1 = M1()*TranslateMatrix(vec3(camera2D.x,camera2D.y,0));
    mat4 m2 = M2()*TranslateMatrix(vec3(camera2D.x,camera2D.y,0));
    gpuProgram.setUniform(m1, "MVP");
    gpuProgram.setUniform(vec3(1,1,1),"color");
    for (int i = 0; i < numberOfEdges1; ++i) {
        glDrawArrays(GL_LINES, i*2, 2);
    }
    gpuProgram.setUniform(m2, "MVP");
    for (int i = 0; i < numberOfEdges2 ; ++i) {
        glDrawArrays(GL_LINES, numberOfEdges1*2+i*2, 2);
    }
    gpuProgram.setUniform(m1, "MVP");
    for (int i = 0; i < numberOfAtoms1; ++i) {
        gpuProgram.setUniform(atoms1[i].color,"color");
        glDrawArrays(GL_TRIANGLE_FAN, (numberOfEdges1+numberOfEdges2)*2+i*CIRCLE_RESOLUTION, (CIRCLE_RESOLUTION-2));
    }
    gpuProgram.setUniform(m2, "MVP");
    for (int i = 0; i < numberOfAtoms2; ++i) {
        gpuProgram.setUniform(atoms2[i].color,"color");
        glDrawArrays(GL_TRIANGLE_FAN, (numberOfEdges1+numberOfEdges2)*2+numberOfAtoms1*CIRCLE_RESOLUTION+i*CIRCLE_RESOLUTION, (CIRCLE_RESOLUTION-2));
    }
    glutSwapBuffers();
}

void onKeyboard(unsigned char key, int pX, int pY) {
    if (key == ' '){
        camera2D = vec2(0,0);
        delete moleculePair;
        moleculePair = new MoleculePair();
    }else if (key == 's'){
        camera2D.x+=0.1;
    }else if(key == 'd'){
        camera2D.x-=0.1;
    }else if(key == 'x'){
        camera2D.y+=0.1;
    }else if(key == 'e'){
        camera2D.y-=0.1;
    }glutPostRedisplay();
}

void onKeyboardUp(unsigned char key, int pX, int pY) { }
void onMouseMotion(int pX, int pY) {}
void onMouse(int button, int state, int pX, int pY) {}

void onIdle() {
    long time = glutGet(GLUT_ELAPSED_TIME);
    long timeBehind = (time - lastUpdate);
    for (int t = 0; t < timeBehind; t+=10){
        moleculePair->Calculate();
        moleculePair->velocity1 = moleculePair->velocity1 + moleculePair->acceleration1 * 0.01;
        moleculePair->velocity2 = moleculePair->velocity2 + moleculePair->acceleration2 * 0.01;
        wTranslate1 = wTranslate1 + moleculePair->velocity1 * 0.01;
        wTranslate2 = wTranslate2 + moleculePair->velocity2 * 0.01;
        moleculePair->angular_vel1 += moleculePair->angular_acc1 * 0.01;
        moleculePair->angular_vel2 += moleculePair->angular_acc2 * 0.01;
        phi1 += moleculePair->angular_vel1 * 0.01f;
        phi2 += moleculePair->angular_vel2 * 0.01f ;
        vec4 tmp1 = vec4(moleculePair->position1.x, moleculePair->position1.y, 0,1) * TranslateMatrix(vec3(wTranslate1.x, wTranslate1.y, 0));
        vec4 tmp2 = vec4(moleculePair->position2.x, moleculePair->position2.y, 0,1) * TranslateMatrix(vec3(wTranslate2.x, wTranslate2.y, 0));
        moleculePair->position1 = vec2(tmp1.x, tmp1.y);
        moleculePair->position2 = vec2(tmp2.x, tmp2.y);
    }
    lastUpdate = time;
    glutPostRedisplay();
}