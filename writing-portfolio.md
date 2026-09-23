---
layout: page
title: Writing Portfolio
permalink: /writing-portfolio/
---

Welcome to my portfolio of written pieces. Below you will find examples of my technical and scientific writing.

<ul>
  {% for post in site.categories.technical-writing %}
    <li>
      <span class="post-meta">{{ post.date | date: "%b %-d, %Y" }}</span>
      <h3><a class="post-link" href="{{ post.url | relative_url }}">{{ post.title | escape }}</a></h3>
    </li>
  {% endfor %}
</ul>
