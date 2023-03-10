package model.strategy;//
//
//  Generated by StarUML(tm) Java Add-In
//
//  @ Project : Untitled
//  @ File Name : strategy.NoLearn.java
//  @ Date : 2022. 03. 23.
//  @ Author : 
//
//


import model.Virologist;
import test.Tester;

/**
 * Tanulási stratégia, amely megakadályozza, hogy egy virológus megtanulhasson genetikai kódot az adott mezőről.
 */
public class NoLearn implements ILearnStr
{
	/**
	 * A stratégia alkalmazásakor hívott metódus.
	 * Nem végez semmilyen műveletet.
	 * @param v A virológus, aki tanulni próbál.
	 */
	public void Learn(Virologist v)
	{
		Tester.methodStart(new Object(){}.getClass().getEnclosingMethod());

		Tester.methodEnd(new Object(){}.getClass().getEnclosingMethod());
	}
}
