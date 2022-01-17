/* 
AP Physics Lab 1A 
Collision Simulation
Sarah Ali, George Paraschov, Arthur Xu
SPH4U0
*/

// define standard mass, radius, time increments, spring constant, and coefficient of friction
#define M 0.017
#define R 0.055
#define DT 0.00001
#define K 426
#define C 0

// libraries
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>

// using standard namespace
using namespace std;

// Vector Class
class Vector
{
public:
    double x, y, m, theta; // variables for components, magnitude, and angle

    Vector() // default constructor for 0 vector
    {
        x = 0;
        y = 0;
        m = 0;
        theta = 0;
    }

    Vector(double a, double b) // vector constructor with user-input rectangular components
    {
        x = a;
        y = b;
        m = sqrt(x * x + y * y); // calculates magnitude
        theta = atan2(y, x);     // calculates angle in radians
    }
};

// overrides output operator to print vectors in readable format (x,y)
ostream &operator<<(ostream &strm, const Vector &v)
{
    return strm << "(" << to_string(v.x) << "," << to_string(v.y) << ")";
}

// Ball class
class Ball
{
public:
    Vector s, v, a;        // linear kinematics vectors
    double m, w, I, alpha; // mass and rotational variables: angular velocity, rotational inertia, angular acceleration

    Ball(Vector s0 = Vector(), Vector v0 = Vector(), double m0 = M, double w0 = 0) //constructor for ball that assigns user input or default values
    {
        s = s0;
        v = v0;
        m = m0;
        w = w0;
    }

    void applyForce(Vector F, Vector compress, double torque) // method that applies a force to effect motion of ball
    {
        a = Vector(F.x / m, F.y / m); //calculate acceleration using N2: a = F/m

        // calculate rotational inertia for an ellipsoid: I = 0.2M(a^2 + b^2) where a and b are semi-minor axis
        I = 0.2 * m * (pow(R - compress.m, 2) + pow(R, 2));
        // calculates angular acceleration using alpha = torque / I
        alpha = torque / I;

        // updates position with kinematics equation s = s0 + v0 * t + 0.5 * a * t^2
        s = Vector(s.x + v.x * DT + 0.5 * a.x * DT * DT, s.y + v.y * DT + 0.5 * a.y * DT * DT);
        // updates linear velocity with kinematics equation v = v0 + a * t
        v = Vector(v.x + a.x * DT, v.y + a.y * DT);
        // updates angular velocity with rotational kinematics equation w = w0 + alpha * t
        w += alpha * DT;
    }
};

// driver function
int main()
{
    double sx, sy, vx, vy, w, m, T; // temporary input values
    double torque;
    Vector spring = Vector();   // spring force vector
    Vector friction = Vector(); // friction force vector
    Vector force = Vector();    // net force on balls

    cout << "Enter the x component of ball 1's position: ";
    cin >> sx;
    cout << "Enter the y component of ball 1's position: ";
    cin >> sy;
    cout << "Enter the x component of ball 1's velocity: ";
    cin >> vx;
    cout << "Enter the y component of ball 1's velocity: ";
    cin >> vy;
    cout << "Enter the angular velocity of ball 1: ";
    cin >> w;
    cout << "Enter the mass of ball 1 (0 for default 0.017 kg): ";
    cin >> m;
    if (m == 0) // allows user to choose default mass
    {
        m = M;
    }

    // defines first object using user input
    Ball b1 = Ball(Vector(sx, sy), Vector(vx, vy), m, w);

    cout << "Enter the x component of ball 2's position: ";
    cin >> sx;
    cout << "Enter the y component of ball 2's position: ";
    cin >> sy;
    cout << "Enter the x component of ball 2's velocity: ";
    cin >> vx;
    cout << "Enter the y component of ball 2's velocity: ";
    cin >> vy;
    cout << "Enter the angular velocity of ball 2: ";
    cin >> w;
    cout << "Enter the mass of ball 2 (0 for default 0.017 kg): ";
    cin >> m;
    if (m == 0) // allows user to choose default mass
    {
        m = M;
    }

    // defines second object using user input
    Ball b2 = Ball(Vector(sx, sy), Vector(vx, vy), m, w);

    cout << "Enter the time limit: ";
    cin >> T;

    // opens output csv file
    string output = "output.csv";
    ofstream f(output);

    // header for output file
    f << "t" << "," << "s1" << "," << "s2" << "," << "v1" << "," << "v2" << "," << "w1" << "," << "w2" << "," << "\n";

    // loop for simulation during time interval with dt increments
    for (double t = 0; t < T; t += DT)
    {
        // outputs values to csv file
        f << to_string(t) << "," << b1.s << "," << b2.s << "," << b1.v << "," << b2.v << "," << b1.w << "," << b2.w << ","
          << "\n";

        // calculates distance vector from ball 1 to ball 2
        Vector distance = Vector(b1.s.x - b2.s.x, b1.s.y - b2.s.y);
        // calculates compression using radii and distance vector
        Vector compression = Vector((2 * R - distance.m) * cos(distance.theta), (2 * R - distance.m) * sin(distance.theta));

        if (compression.m > 0) // computes motion during collision
        {
            // calculates spring force using F = kx
            spring = Vector(compression.x * K, compression.y * K);
            // friction force is calculated as friction = coeff of friction * normal force where the spring force is the normal force
            friction = Vector(-spring.y * C, spring.x * C);
            // torque is taken as the 2D pseudo-cross product: r x F
            // r is the radius to the point of contact and F is only the force of friction since the spring force is perpendicular
            torque = (R * cos(compression.theta) - compression.m / 2.0) * friction.y - (R * sin(compression.theta) - compression.m / 2.0) * friction.x;

            // checks if the spring force would be parallel to the direction of friction since the spring
            // force is the normal force and friction is antiparallel to the velocity
            if (distance.y / distance.x - b1.v.y / b1.v.x < pow(1, -5))
            {
                // if they are parallel friction does not affect the translational motion
                // the translational force is only the spring force
                force = spring;
            }
            else
            {
                // the net translational force is the spring force and friction force
                force = Vector(spring.x + friction.x, spring.y + friction.y);
            }
        }
        else
        {
            // resets values to 0 if balls are not in contact
            force = Vector();
            compression = Vector();
            torque = 0;
        }

        // updates kinematic values by applying spring force to proceed to next step in simulation
        b1.applyForce(force, compression, torque);
        b2.applyForce(Vector(-force.x, -force.y), compression, -torque);
    }

    // closes output file
    f.close();
}