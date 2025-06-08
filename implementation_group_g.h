/* Implementation - Group G: Pretty Printing Functions */

#ifndef FIXSTRING_NO_STDIO

/*---------------------------------------------------------------------------*/
/* PRETTY PRINTING FUNCTIONS                                                 */
/*---------------------------------------------------------------------------*/

/* Base formatting function for fix16_t values */
FIXMATH_FUNC_ATTRS void print_fix16_t(FILE *stream, fix16_t value, uint_fast8_t width, uint_fast8_t decimals)
{
    char buf[13];
    fix16_to_str(value, buf, decimals);
    
    uint_fast8_t len = strlen(buf);
    if (len < width)
    {
        width -= len;
        while (width-- > 0)
            fputc(' ', stream);
    }

    fputs(buf, stream);
}

/* Pretty print matrix */
FIXMATH_FUNC_ATTRS void print_mf16(FILE *stream, const mf16 *matrix)
{
    if (matrix->errors)
    {
        fprintf(stream, "MATRIX ERRORS: %d\n", matrix->errors);
    }
    
    int row, column;
    for (row = 0; row < matrix->rows; row++)
    {
        for (column = 0; column < matrix->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            print_fix16_t(stream, value, 9, 4);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
    }
}

#endif /* FIXSTRING_NO_STDIO */ 