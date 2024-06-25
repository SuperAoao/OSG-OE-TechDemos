uniform int window_width_half;
uniform int window_height_half;

varying float window_pos_x; // doesn't allow varying int
varying float window_pos_y;

void main()
{
    gl_Position = ftransform();
    window_pos_x = window_width_half * gl_Position.x / gl_Position.w + window_width_half;
    window_pos_y = window_height_half * gl_Position.y / gl_Position.w + window_height_half;
}
