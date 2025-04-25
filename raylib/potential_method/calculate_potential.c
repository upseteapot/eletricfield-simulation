#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#define THREADS 4


typedef struct
{
    float value;
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

void AddPointPotentialArray(potential_array_t *array, float value, int x, int y)
{
    if (array->size == array->capacity)
        return;
    
    for (size_t i=0; i < array->size; i++)
        if (array->data[i].x == x && array->data[i].y == y)
            return;
    
    array->data[array->size++] = (potential_t){ value, x, y };
}

void AddLinePotentialArray(potential_array_t *array, float value, int x0, int y0, int x1, int y1)
{
    float m = (float)(y1 - y0)/(x1 - x0);
    for (int x=x0; x <= x1; x++)
        AddPointPotentialArray(array, value, x, y0 + m * (x - x0));
}

void AssignPotentialToEscalarField(potential_array_t *array, float *escalar_field, int columns)
{
    for (size_t i=0; i < array->size; i++)
    {
        float value = array->data[i].value;
        int x = array->data[i].x;
        int y = array->data[i].y;
        escalar_field[x + y * columns] = value;
    }
}

typedef struct
{
    float *source;
    float *dest;
    int y_start;
    int y_stop;
    int rows;
    int columns;
} cps_param_t;

void *CalculatePotentialStep(void *param)
{
    cps_param_t *cps_param = (cps_param_t*)param;
    float *source = cps_param->source;
    float *dest = cps_param->dest;
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
    const float y1 = 0.0014;
    const float y2 = 0.0064;
    const float x1 = 0.0080;
    const float x2 = 0.0220;
    const float screen_distance = 0.225;
    const float margin_y = 0.01;
    const float margin_x = 0.005;

    const float division = 0.0001 / 3;
    const float world_height = (y2 + 2*margin_y) / division;
    const float world_width = (x1 + x2 + screen_distance + 2*margin_x) / division;
    const size_t rows = world_height / division;
    const size_t columns = world_width / division;
    
    const size_t simulation_steps = 10000;
    
    printf("[DEBUG]\n");
    printf("[rows] %zu, [columns] %zu\n", rows, columns);
    
    float eletric_field[2 * rows * columns]; 
    memset(eletric_field, '\0', sizeof(float) * 2 * rows * columns);
    
    //float new_potential[rows * columns];
    //float potential[rows * columns];
    //memset(new_potential, '\0', sizeof(float) * rows * columns); 
    //memset(potential, '\0', sizeof(float) * rows * columns); 
    //
    //potential_array_t fixed_potential;
    //CreatePotentialArray(&fixed_potential, rows * columns);

    //AddLinePotentialArray(&fixed_potential, 1,
    //    (margin_x)/division,
    //    rows/2 - (y1/2.0)/division,
    //    (margin_x/2.0f + x1)/division,
    //    rows/2 - (y1/2.0)/division
    //);
    //AddLinePotentialArray(&fixed_potential, 1,
    //    (margin_x + x1)/division,
    //    rows/2 - (y1/2.0)/division,
    //    (margin_x/2.0f + x1 + x2)/division,
    //    rows/2 - (y2/2.0)/division
    //);
    //AddLinePotentialArray(&fixed_potential, 0,
    //    (margin_x)/division,
    //    rows/2 + (y1/2.0)/division,
    //    (margin_x + x1)/division,
    //    rows/2 + (y1/2.0)/division
    //);
    //AddLinePotentialArray(&fixed_potential, 0,
    //    (margin_x + x1)/division,
    //    rows/2 + (y1/2.0)/division,
    //    (margin_x + x1 + x2)/division,
    //    rows/2 + (y2/2.0)/division
    //);
    //
    //pthread_t thread_pool[THREADS];
    //cps_param_t parameters[THREADS];
    //
    //printf("Evaluating...\n");
    //
    //for (size_t s=0; s < simulation_steps; s++)
    //{
    //    AssignPotentialToEscalarField(&fixed_potential, potential, columns);

    //    for (size_t t=0; t < THREADS; t++)
    //    {
    //        parameters[t].source = potential;
    //        parameters[t].dest = new_potential;
    //        parameters[t].y_start = ((float)(t) / THREADS) * rows;
    //        parameters[t].y_stop = ((float)(t + 1) / THREADS) * rows;
    //        parameters[t].rows = rows;
    //        parameters[t].columns = columns;
    //        pthread_create(thread_pool+t, NULL, CalculatePotentialStep, (void *)(parameters+t));
    //    }

    //    for (int t=0; t < THREADS; t++) 
    //        pthread_join(thread_pool[t], NULL);

    //    memcpy(potential, new_potential, sizeof(float) * rows * columns); 
    //}

    //AssignPotentialToEscalarField(&fixed_potential, potential, columns);
    //
    //for (int y=0; y < rows; y++)
    //for (int x=0; x < columns; x++)
    //{
    //    int n_y = (y == 0) ? (rows - 1) : (y - 1);
    //    int n_x = (x == 0) ? (columns - 1) : (x - 1);

    //    eletric_field[2*(x + y * columns)] = (potential[n_x + y * columns] - potential[x + y * columns]) / (2*division);
    //    eletric_field[2*(x + y * columns) + 1] = (potential[x + n_y * columns] - potential[x + y * columns]) / (2*division);
    //}
    //
    //printf("Saving potential...\n");
    //
    //FILE *field_file = fopen("normalized_field", "w");
    //FILE *potential_file = fopen("normalized_potential", "w");
    //
    //fwrite((void *)&rows, sizeof(size_t), 1, field_file);
    //fwrite((void *)&columns, sizeof(size_t), 1, field_file);
    //fwrite((void *)&rows, sizeof(size_t), 1, potential_file);
    //fwrite((void *)&columns, sizeof(size_t), 1, potential_file);

    //fwrite((void *)&eletric_field, sizeof(float) * 2, rows * columns, field_file);
    //fwrite((void *)&potential, sizeof(float) * 2, rows * columns, field_file);

    //fclose(field_file);
    //fclose(potential_file);
    //
    //printf("Finished!\n");

    return 0;
}

