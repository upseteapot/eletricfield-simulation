#include <stdio.h>
#include <raylib.h>
#include <stdbool.h>
#include <pthread.h>
#include <string.h>
#include <math.h>

#define THREADS 4

static const float scale = 0.01;
static const float e0 = 8.85 * pow(10, -12);

void interpolate_color(Color *color, float x, float min, float max, Color start, Color end)
{
    float a = (x - min) / (max - min);
    color->r = a * (end.r - start.r) + start.r;
    color->g = a * (end.g - start.g) + start.g;
    color->b = a * (end.b - start.b) + start.b;
    color->a = a * (end.a - start.a) + start.a;
}

typedef struct
{
    int x;
    int y;
    float value;
} charge_t;

void add_charge(charge_t *charges, size_t *size, int x, int y, float value)
{
    for (size_t i=0; i < *size; i++)
    {
        if (charges[i].x == x && charges[i].y == y)
        {
            charges[i].value = value;
            return;
        }
    }

    charges[(*size)] = (charge_t){ x, y, value };
    (*size) += 1;
}

typedef struct
{
    float *field;
    charge_t *charges;
    size_t size;
    int y0;
    int y1;
    int rows;
    int columns;
} simulate_parameter_t;

void *simulate(void *param)
{
    simulate_parameter_t *simulate_param = (simulate_parameter_t *)param; 

    float *field = simulate_param->field;
    charge_t *charges = simulate_param->charges;
    size_t size = simulate_param->size;
    int y0 = simulate_param->y0;
    int y1 = simulate_param->y1;
    int rows = simulate_param->rows;
    int columns = simulate_param->columns;
    
    printf("y0:%d, y1:%d, cols:%d, rows:%d\n", y0, y1, columns, rows);

    for (int y=y0; y < y1; y++)
    for (int x=0; x < columns; x++)
    {
        for (size_t i=0; i < size; i++)
            if (!(y == charges[i].y && x == charges[i].x))
                field[x + y*columns] += 
                    1 / (4 * PI * e0) * (charges[i].value) 
                    / (scale * sqrt(pow(y - charges[i].y, 2) + pow(x - charges[i].x, 2)));
    }
}

int main(void)
{
    // TODO: Better circle algorithm.
    // TODO: two methods: calculate potencial (average) then calculate the eletric field (div).
    //                    calculate eletric field by adding all dQ contribution.
    // TODO: Draw with shaders
    InitWindow(0,0,"Eletric Field Simulation");
    
    const int width = GetScreenWidth();
    const int height = GetScreenHeight();
    
    const int resolution = 3;
    const int columns = width / resolution;
    const int rows = height / resolution;
    const int x_offset = (width - columns * resolution) / 2;
    const int y_offset = (height - rows * resolution) / 2;
    
    float field[rows * columns];
    charge_t charges[rows * columns];
    size_t size = 0;
    memset(field, '\0', sizeof(float) * rows * columns);
    memset(charges, '\0', sizeof(charge_t) * rows * columns);
    
    float max_field = -INFINITY;
    float min_field = INFINITY;
    
    bool draw_grid = false;
    bool run_simulation = false;
    bool draw_field = false;

    Color cell_color;
    Color min_color = {0, 0, 255, 255};
    Color max_color = {255, 0, 0, 255};
    
    Vector2 center;
    bool center_defined = false;
    
    while (!WindowShouldClose())
    {
        Vector2 window_mouse_pos = GetMousePosition();
        Vector2 mouse_pos = { (window_mouse_pos.x - x_offset) / resolution, (window_mouse_pos.y - y_offset) / resolution };
        
        BeginDrawing();
        ClearBackground(BLACK);
        
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && !center_defined)
        {
            center = mouse_pos;
            center_defined = true;
        }
        
        else if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && center_defined)
        {
            float radius = sqrt(pow(center.x - mouse_pos.x, 2) + pow(center.y - mouse_pos.y, 2));
            float dtheta = 0.00001f;
            
            printf("%f, %f\n", radius, dtheta);

            for (float theta=0; theta < 2 * PI; theta += dtheta)
            {
                int x1 = center.x + radius * cos(theta);
                int y1 = center.y + radius * sin(theta);
                add_charge(charges, &size, x1, y1, 1);
            }
            
            center_defined = false;
        }
        
        if (IsMouseButtonPressed(MOUSE_RIGHT_BUTTON))
        {
            int x = (GetMouseX() - x_offset) / resolution;
            int y = (GetMouseY() - y_offset) / resolution;

            printf("[%d, %d]\n", x, y);
            printf("\t\t%f\n", field[x + (y-1)*columns]);
            printf("%f\t%f\t%f\n", field[x-1 + y*columns], field[x + y*columns], field[x+1 + y*columns]);
            printf("\t\t%f\n", field[x + (y+1)*columns]);
        }
        
        run_simulation = IsKeyPressed(KEY_SPACE);

        // Draw field.
        if (draw_field) 
        {
            for (int y=0; y < rows; y++)
                for (int x=0; x < columns; x++)
                {
                    interpolate_color(&cell_color, field[x + y*columns], min_field, max_field, min_color, max_color);
                    DrawRectangle(x_offset + x*resolution, y_offset + y*resolution, resolution, resolution, cell_color);
                }
        }
        
        // Draw charges.
        for (size_t i=0; i < size; i++)
        {
            DrawRectangle(x_offset + charges[i].x*resolution, y_offset + charges[i].y*resolution, resolution, resolution, WHITE);
        }
        
        // Draw grid.
        if (draw_grid)
        {
            for (int x=0; x <= columns; x++)
                DrawLine(x_offset + x*resolution, y_offset, x_offset + x*resolution, y_offset + rows*resolution, GRAY);
            for (int y=0; y <= rows; y++)
                DrawLine(x_offset, y_offset + y*resolution, x_offset + columns*resolution, y_offset + y*resolution, GRAY);
        }

        if (run_simulation)
        {
            draw_field = true;

            memset(field, '\0', sizeof(float) * rows * columns);

            pthread_t threads[THREADS];
            simulate_parameter_t parameters[THREADS];
            
            for (int i=0; i < THREADS; i++) 
            {
                parameters[i].field = field;
                parameters[i].charges = charges;
                parameters[i].size = size;
                parameters[i].y0 = ((float)(i) / THREADS) * rows;
                parameters[i].y1 = ((float)(i + 1) / THREADS) * rows;
                parameters[i].rows = rows;
                parameters[i].columns = columns;
                pthread_create(threads+i, NULL, simulate, (void *)(parameters+i));
            }
            
            for (int i=0; i < THREADS; i++) 
            {
                pthread_join(threads[i], NULL);
            }

            for (int y=0; y < rows; y++)
                for (int x=0; x < columns; x++)
                {
                    max_field = fmax(max_field, field[x + y*columns]);
                    min_field = fmin(min_field, field[x + y*columns]);
                }

            printf("max:%f\n", min_field);
            printf("min:%f\n", max_field);

            run_simulation = false;
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}

