#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <raylib.h>
#define THREADS 4


// Interpolate two colors. 
void interpolate_color(Color *color, double x, double min, double max, Color start, Color end)
{
    double a = (x - min) / (max - min);
    color->r = a * (end.r - start.r) + start.r;
    color->g = a * (end.g - start.g) + start.g;
    color->b = a * (end.b - start.b) + start.b;
    color->a = a * (end.a - start.a) + start.a;
} 

// Potential array.
typedef struct
{
    double value;
    int x;
    int y;
} potential_t;

typedef struct
{
    potential_t *data;
    size_t size;
    size_t capacity;
} potential_array_t;

void CreatePotentialArray(potential_array_t *array, size_t capacity)
{
    array->data = (potential_t *)malloc(sizeof(potential_t) * capacity);
    array->size = 0;
    array->capacity = capacity;
}

void AddPointPotentialArray(potential_array_t *array, double value, int x, int y)
{
    if (array->size == array->capacity)
        return;
    
    for (size_t i=0; i < array->size; i++)
        if (array->data[i].x == x && array->data[i].y == y)
            return;
    
    array->data[array->size++] = (potential_t){ value, x, y };
}

void AddLinePotentialArray(potential_array_t *array, double value, int x0, int y0, int x1, int y1)
{
    double m = (double)(y1 - y0)/(x1 - x0);
    for (int x=x0; x <= x1; x++)
        AddPointPotentialArray(array, value, x, y0 + m * (x - x0));
}

void AssignPotentialToEscalarField(potential_array_t *array, double *escalar_field, int columns)
{
    for (size_t i=0; i < array->size; i++)
    {
        double value = array->data[i].value;
        int x = array->data[i].x;
        int y = array->data[i].y;
        escalar_field[x + y * columns] = value;
    }
}

// Potential calculation step.
typedef struct
{
    double *source;
    double *dest;
    int y_start;
    int y_stop;
    int rows;
    int columns;
} cps_param_t;

void *CalculatePotentialStep(void *param)
{
    cps_param_t *cps_param = (cps_param_t*)param;
    double *source = cps_param->source;
    double *dest = cps_param->dest;
    int y_start = cps_param->y_start;
    int y_stop = cps_param->y_stop;
    int rows = cps_param->rows;
    int columns = cps_param->columns;
    
    int x1, x2, y1, y2;
    
    for (int y=y_start; y < y_stop; y++)
    for (int x=0; x < columns; x++)
    {
        // Calculate new potential.
        y1 = (y == 0) ? (rows - 1) : (y - 1);
        y2 = (y == rows - 1) ? 0 : (y + 1);
        x1 = (x == 0) ? (columns - 1) : (x - 1);
        x2 = (x == columns - 1) ? 0 : (x + 1);
        
        dest[x + y * columns] = (1. / 4.) * 
            (source[x1 + y  * columns] +
             source[x  + y1 * columns] +
             source[x2 + y  * columns] +
             source[x  + y2 * columns]);    
    }
}


int main(void)
{
    // SIMULATION/WORLD SETUP
    // Physical constants. 
    const double electron_charge = -1.602176565 * pow(10, -19);
    const double electron_mass = 9.1093837 * pow(10, -31);
    
    // Apparatus dimensions in meters.
    const double y1 = 0.0014;
    const double y2 = 0.0064;
    const double x1 = 0.0080;
    const double x2 = 0.0220;
    const double screen_distance = 0.225;

    const double unity_size = 0.0001; // Correspond to the size of 1 world unity in meters.
    const double scale = 1/unity_size; // Scale factor (meters to world unities)
    const double division = 2; // Correspond to the "dV" size in world unity.
    const double world_height = (y2 + 0.03) * scale; // World height in world unities.
    const double world_width = (x1 + x2 + screen_distance + 0.01) * scale; // World width in world unities.
    
    // Grid size.
    const int rows = world_height / division;
    const int columns = world_width / division;
    
    // Eletric field vector field and magnitude array.
    Vector2 eletric_field[rows * columns]; 
    double  eletric_field_magnitude[rows * columns];
    double  eletric_field_max;
    double  eletric_field_min;
    memset(eletric_field, '\0', sizeof(Vector2) * rows * columns);
    memset(eletric_field_magnitude, '\0', sizeof(double) * rows * columns);
    
    // New potential and potential escalar field.
    double new_potential[rows * columns];
    double potential[rows * columns];
    double potential_max;
    double potential_min;
    memset(new_potential, '\0', sizeof(double) * rows * columns); 
    memset(potential, '\0', sizeof(double) * rows * columns); 
    
    potential_array_t fixed_potential;
    CreatePotentialArray(&fixed_potential, rows * columns);
    
    // Thread for parallel processing the potential field calculation.
    pthread_t thread_pool[THREADS];
    cps_param_t parameters[THREADS];
    
    // Create apparatus in world.
    const float potential_difference = 14;
    const double V_acc = 752;

    AddLinePotentialArray(&fixed_potential, potential_difference,
        0,
        (world_height + (-y1/2 - y2/2)*(scale))/division,
        x1*(scale/division),
        (world_height + (-y1/2 - y2/2)*(scale))/division
    );
    AddLinePotentialArray(&fixed_potential, potential_difference,
        x1*(scale/division),
        (world_height + (-y1/2 - y2/2)*(scale))/division,
        (x1 + x2)*(scale/division),
        (world_height + (-y2)*(scale))/division
    );
    AddLinePotentialArray(&fixed_potential, 0,
        0,
        (world_height + (y1/2 - y2/2)*(scale))/division,
        x1*(scale/division),
        (world_height + (y1/2 - y2/2)*(scale))/division
    );
    AddLinePotentialArray(&fixed_potential, 0,
        x1*(scale/division),
        (world_height + (y1/2 - y2/2)*(scale))/division,
        (x1 + x2)*(scale/division),
        (world_height/division)
    );
    
    // Create particle's initial state.
    Vector2 particle_pos = (Vector2){ 
        0.0f,
        world_height - y2*scale/2.0f
    };
    
    Vector2 particle_vel = (Vector2){ 
        sqrt(fabs(2 * electron_charge * V_acc / electron_mass))*scale, 
        0.0f 
    };
    
    // Simulation configuration.
    const double dt = pow(10, -13);
    const size_t particle_simulation_steps = 1000;
    const size_t simulation_steps = 5000;
    
    // GRAPHICS/WINDOW SETUP
    InitWindow(0, 0,"Eletric Field Simulation");
    const int window_width = GetScreenWidth();
    const int window_height = GetScreenHeight();
    const double scale_to_pixels = window_width / world_width; // World unities to pixels.
    
    bool simulate_potential = true;
    bool simulate_particle = false;
    bool draw_field = false;
    bool draw_potential = false;
    bool draw_fixed_potential = true;
    bool draw_particle = true;
    
    Color cell_color;
    Color potential_max_color = {255, 0, 0, 255};
    Color potential_min_color = {0, 0, 255, 255};
    Color field_max_color = {0, 255, 0, 255};
    Color field_min_color = {0, 0, 255, 255};
    
    // DEBUG
    printf("\n[world_width] %f, [world_height] %f\n", world_width, world_height);
    printf("[scale] %f, [scale_to_pixels] %f\n", scale, scale_to_pixels);
    printf("[particle_x] %f, [particle_y] %f\n", particle_pos.x, particle_pos.y);
    printf("[s_particle_x] %f, [s_particle_y] %f\n", scale_to_pixels * particle_pos.x, scale_to_pixels * particle_pos.y);
    printf("[fixed_potential_size] %zu\n", fixed_potential.size);
    printf("[V_acc] %f, [dp] %f\n", V_acc, potential_difference);
    
    // START SIMULATION
    while (!WindowShouldClose())
    {
        if (IsKeyPressed(KEY_SPACE))
            simulate_potential = !simulate_potential;
        
        if (IsKeyPressed(KEY_P))
            simulate_particle = !simulate_particle;
        
        if (IsKeyPressed(KEY_S))
        {
            draw_field = !draw_field;
            draw_potential = !draw_field;
        }
        
        // Simulate potential.
        if (simulate_potential)
        {
            // Calculate potential.
            for (size_t s=0; s < simulation_steps; s++)
            {
                AssignPotentialToEscalarField(&fixed_potential, potential, columns);
                
                for (size_t t=0; t < THREADS; t++)
                {
                    parameters[t].source = potential;
                    parameters[t].dest = new_potential;
                    parameters[t].y_start = ((double)(t) / THREADS) * rows;
                    parameters[t].y_stop = ((double)(t + 1) / THREADS) * rows;
                    parameters[t].rows = rows;
                    parameters[t].columns = columns;
                    pthread_create(thread_pool+t, NULL, CalculatePotentialStep, (void *)(parameters+t));
                }

                for (int t=0; t < THREADS; t++) 
                    pthread_join(thread_pool[t], NULL);

                memcpy(potential, new_potential, sizeof(double) * rows * columns); 
            }
            
            AssignPotentialToEscalarField(&fixed_potential, potential, columns);
            
            // Calculate eletric field.
            for (int y=0; y < rows; y++)
            for (int x=0; x < columns; x++)
            {
                int n_y = (y == 0) ? (rows - 1) : (y - 1);
                int n_x = (x == 0) ? (columns - 1) : (x - 1);

                eletric_field[x + y * columns] = (Vector2){
                    (potential[n_x + y * columns] - potential[x + y * columns]) / (2*unity_size),
                    (potential[x + n_y * columns] - potential[x + y * columns]) / (2*unity_size)
                }; 

                eletric_field_magnitude[x + y * columns] = sqrt(
                    pow(eletric_field[x + y * columns].x, 2) + 
                    pow(eletric_field[x + y * columns].y, 2)
                );
            }
            
            FILE *field_file = fopen("field", "w");
            FILE *potential_file = fopen("potential", "w");
            int x_start = 0;
            int x_end = (x1 + x2)*(scale/division);
            int y = (world_height - y2*scale/2.0)/division;
            for (int x=x_start; x <= x_end; x++)
            {
                fprintf(field_file, "%f %f\n", 
                    x * division * unity_size, 
                    eletric_field_magnitude[x + y * columns]
                );
                
                fprintf(potential_file, "%f %f\n", 
                    x * division * unity_size, 
                    potential[x + y * columns]
                );
            
            }
            fclose(field_file);
            fclose(potential_file);
            
            potential_max = -INFINITY;
            potential_min =  INFINITY;
            eletric_field_max = -INFINITY;
            eletric_field_min =  INFINITY;
            double value;

            for (int y=0; y < rows; y++)
            for (int x=0; x < columns; x++)
            {
                value = potential[x + y * columns];
                potential_max = fmax(value, potential_max);
                potential_min = fmin(value, potential_min);
                
                value = eletric_field_magnitude[x + y * columns];
                eletric_field_max = fmax(value, eletric_field_max);
                eletric_field_min = fmin(value, eletric_field_min);
            }
            
            printf("[max_field] %f; [min_field] %f\n", eletric_field_max, eletric_field_min);
            printf("[potential_max] %f; [potential_min] %f\n\n", potential_max, potential_min);
            
            simulate_potential = false;
            draw_potential = true;
            simulate_particle = true;
        }
        
        if (simulate_particle)
        {
            for (size_t i=0; i < particle_simulation_steps; i++)
            {
                int x = particle_pos.x / division; 
                int y = particle_pos.y / division; 

                Vector2 field;
                if (x < columns && y < rows)
                    field = eletric_field[x + y * columns];
                else
                    field = (Vector2){ 0, 0 };

                particle_vel.x += ((field.x * electron_charge / electron_mass) * dt) * scale;
                particle_vel.y += ((field.y * electron_charge / electron_mass) * dt) * scale;

                particle_pos.x += particle_vel.x * dt;
                particle_pos.y += particle_vel.y * dt;
                
                //printf("[field] %f %f\n", field.x, field.y);
                //printf("[pos] %f %f\n", particle_pos.x, particle_pos.y);
                //printf("[dpos] %f %f\n", particle_vel.x * dt, particle_vel.y * dt);
                //printf("[vel] %f %f\n\n", particle_vel.x, particle_vel.y);

                if (particle_pos.x >= (x1 + x2 + screen_distance)*scale)
                {
                    printf("[pos] %f %f\n", particle_pos.x, particle_pos.y);
                    printf("[dy] %f\n", (world_height - y2*scale/2.0 - particle_pos.y) * unity_size);
                    simulate_particle = false;
                    break;
                }
            }
        }
        
        BeginDrawing();
        ClearBackground(BLACK);
        
        // Draw eletric field.
        if (draw_field)
        {
            for (int y=0; y < rows; y++)
            for (int x=0; x < columns; x++)
            {
                interpolate_color(
                    &cell_color, 
                    eletric_field_magnitude[x + y * columns],
                    eletric_field_min,
                    eletric_field_max, 
                    field_min_color, 
                    field_max_color
                );
                
                DrawRectangle(
                    x * division * scale_to_pixels, 
                    y * division * scale_to_pixels, 
                    division * scale_to_pixels, 
                    division * scale_to_pixels,
                    cell_color
                );
            }
        }
        
        // Draw potential
        if (draw_potential)
        {
            for (int y=0; y < rows; y++)
            for (int x=0; x < columns; x++)
            {
                interpolate_color(
                    &cell_color, 
                    potential[x + y * columns],
                    potential_min,
                    potential_max, 
                    potential_min_color, 
                    potential_max_color
                );
                
                DrawRectangle(
                    x * division * scale_to_pixels, 
                    y * division * scale_to_pixels, 
                    division * scale_to_pixels, 
                    division * scale_to_pixels,
                    cell_color
                );
            }
        }
        
        // Draw fixed potential
        if (draw_fixed_potential)
        {
            for (size_t i=0; i < fixed_potential.size; i++)
            {
                int x = fixed_potential.data[i].x;
                int y = fixed_potential.data[i].y;
                
                DrawRectangle(
                    x * division * scale_to_pixels, 
                    y * division * scale_to_pixels, 
                    division * scale_to_pixels, 
                    division * scale_to_pixels,
                    WHITE
                );
            }
        }

        // Draw particle
        if (draw_particle)
        {
            DrawCircle(
                particle_pos.x * scale_to_pixels,
                particle_pos.y * scale_to_pixels, 
                2, 
                WHITE
            );
        }
        
        DrawLine(
            0,
            (world_height - y2*scale/2.0) * scale_to_pixels,
            (x1 + x2 + screen_distance) * scale * scale_to_pixels,
            (world_height - y2*scale/2.0) * scale_to_pixels,
            GREEN
        );

        EndDrawing(); 
    }

    printf("\n");
    CloseWindow();

    return 0;
}

