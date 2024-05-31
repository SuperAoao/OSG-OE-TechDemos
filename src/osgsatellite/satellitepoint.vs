uniform int window_width_half;
uniform int window_height_half;

varying int window_pos_x;
varying int window_pos_y;

void main()
{
    gl_Position = ftransform();
    window_pos_x = int(window_width_half * gl_Position.x / gl_Position.w + window_width_half);
    window_pos_y = int(window_height_half * gl_Position.y / gl_Position.w + window_height_half);
}
