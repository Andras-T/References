package Tetris.Operation;

/**
 * 'T' alakú tetramino
 */
public class ShapeT extends Tetramino{

    /**
     * Létrehoz egy 'T' alakú tetramino-t
     */
    public ShapeT(){
        super();
        int[][] shape = new int[][]{
                {0, 1, 0, 0},
                {0, 1, 1, 0},
                {0, 1, 0, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeT copy(){
        return new ShapeT();
    }

}
