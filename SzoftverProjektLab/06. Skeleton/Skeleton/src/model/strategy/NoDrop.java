package model.strategy;//
//
//  Generated by StarUML(tm) Java Add-In
//
//  @ Project : Untitled
//  @ File Name : strategy.NoDrop.java
//  @ Date : 2022. 03. 23.
//  @ Author : 
//
//


import model.Virologist;
import model.equipments.Equipment;
import test.Tester;

/**
 * Eldobási stratégia, amely megakadályozza, hogy a virológus eldobjon egy felszerelést.
 */
public class NoDrop implements IDropStr
{
	/**
	 * A stratégia alkalmazásakor hívott metódus.
	 * Nem végez semmilyen műveletet.
	 * @param v A virológus, aki el próbál dobni egy felszerelést.
	 * @param e A felszerelés, amit el próbál dobni.
	 */
	public void Drop(Virologist v, Equipment e)
	{
		Tester.methodStart(new Object(){}.getClass().getEnclosingMethod());

		Tester.methodEnd(new Object(){}.getClass().getEnclosingMethod());
	}
}
