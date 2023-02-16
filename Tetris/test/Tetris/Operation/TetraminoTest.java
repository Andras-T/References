package Tetris.Operation;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Tetramino-t tesztel
 */
@RunWith(Parameterized.class)
public class TetraminoTest {

    private final Tetramino shape;

    public TetraminoTest(Tetramino alakzat){
        shape = alakzat;
    }

    /**
     * Szín beállítása
     */
    @Test
    public void colorTest(){
        Assert.assertEquals(shape.getColor(), Color.WHITE);
        shape.setColor(Color.BLUE);
        Assert.assertEquals(shape.getColor(),Color.BLUE);
    }

    /**
     * Másolás teszt
     */
    @Test
    public void copyTest(){
        Assert.assertNotEquals(shape, shape.copy());
    }

    /**
     * Esés teszt
     */
    @Test
    public void dropTest(){
        int y1 = shape.getPosition().y;
        shape.dropTetramino();
        Assert.assertNotEquals(y1, shape.getPosition().y);
    }

    /**
     * Paraméterek beállítása
     * @return Paraméterek
     */
    @Parameters
    public static Collection<Object> parameters() {
        ArrayList<Object> params = new ArrayList<>();
        params.add(new ShapeT());
        params.add(new ShapeI());
        params.add(new ShapeJ());
        params.add(new ShapeL());
        params.add(new ShapeO());
        params.add(new ShapeS());
        params.add(new ShapeZ());
        return params;
    }

}