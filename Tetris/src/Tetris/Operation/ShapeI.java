package Tetris.Operation;

/**
 * 'I' alakú tetramino
 */
public class ShapeI extends Tetramino{

    /**
     * Létrehoz egy 'I' alakú tetramino-t
     */
    public ShapeI(){
        super();
        int[][] shape = new int[][]{
                {0, 1, 0, 0},
                {0, 1, 0, 0},
                {0, 1, 0, 0},
                {0, 1, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeI copy(){
        return new ShapeI();
    }

}
