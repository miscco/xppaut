#ifndef XPPAUT_ARRAY_PRINT_H
#define XPPAUT_ARRAY_PRINT_H

/* --- Macros --- */
#define GREYSCALE -1
#define REDBLUE  0
#define ROYGBIV  1

/* --- Types --- */
typedef struct {
  float xmin,xmax,ymin,ymax;
  float xscale,yscale,xoff,yoff;
  float tx,ty,angle,slant;  /* text attributes   */
  float linecol,letx,lety;
  int linewid;
  } DEVSCALE;

/* --- Functions --- */
int array_print(char *filename, char *xtitle, char *ytitle, char *bottom, int nacross, int ndown, int col0, int row0, int nskip, int ncskip, int maxrow, int maxcol, float **data, double zmin, double zmax, double tlo, double thi, int type);
void ps_replot(float **z, int col0, int row0, int nskip, int ncskip, int maxrow, int maxcol, int nacross, int ndown, double zmin, double zmax, int type);
void ps_begin(double xlo, double ylo, double xhi, double yhi, float sx, float sy);
void ps_convert(double x, double y, float *xs, float *ys);
void ps_col_scale(double y0, double x0, double dy, double dx, int n, double zlo, double zhi, int type, float mx);
void ps_boxit(double tlo, double thi, double jlo, double jhi, double zlo, double zhi, char *sx, char *sy, char *sb, int type);
void ps_close(void);
void ps_setline(double fill, int thick);
void ps_put_char(char ch, float *x, float *y);
void ps_text2(char *str, double xr, double yr, int icent);
void ps_line2(double x1r, double y1r, double x2r, double y2r);
void ps_set_text(double angle, double slant, double x_size, double y_size);
void ps_rect(double x, double y, double wid, double len);
void ps_bar(double x, double y, double wid, double len, double fill, int flag);
void ps_rgb_bar(double x, double y, double wid, double len, double fill, int flag, int rgb);
void ps_hsb_bar(double x, double y, double wid, double len, double fill, int flag);

#endif /* XPPAUT_ARRAY_PRINT_H */
